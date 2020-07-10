// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
template<typename V, typename E>
class Lprop {
public:
	using GraphT = CSRGraph<V, E>;
	Lprop(GraphT* G):G(G){}
    double run(int niter);
    int get_count();
    void errorCheck(std::string message);
    void experiment(const int niter,const bool opt_load, const bool onestep, const bool if_inverse, std::string filename);
    int binary_search(E left,E right,std::vector<V>  offsets,std::vector<V> vertices,V target);

    GraphT* G;
protected:
    void preprocess();
    void postprocess();
    double perform_lp(V n,E m,int niter);
    void init_gmem(V n,E m);
    void free_gmem();
    // Attributes
    V* labels;
    int *d_vertices;     // n
    int *d_neighbors;    // m
    E *d_offsets;      // n + 1
    int *d_labels;       // n
    int *d_sec_labels;
    int *d_counter;      // 1

    int *d_adj_labels; //m
    int *d_sec_adj_labels;//m
    bool *d_if_update; //n

    // for huge nodes
    int *d_block_num;
	int* d_block_v;
	int* d_block_id;
	//for medium nodes
	int *d_warp_num;
	int* d_warp_v;
	int* d_warp_id;
	GlobalHT d_tables; //big and huge edges *2
	//for glolbal reduce
	unsigned long long  *d_max_count_huge;
    // for y-x sort
    int *d_neighborkey_out; //m
    int *d_neighborval_out; //m

    struct configure{
        	bool opt_load = false;
    		// whether one step look up, set opt_load = true before set onestep false.
    		bool onestep = true;
    		// whether inverse when load label.
    		bool if_inverse = false;
    		// when to use Boolean array.
    		int changeiter = 8;
        } myconf;
};

template<typename V, typename E>
void Lprop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);
}
template<typename V, typename E>
void Lprop<V, E>::postprocess()
{
    this->free_gmem();
}

template<class V, class E>
void Lprop<V, E>::errorCheck(std::string message){
	auto err = cudaGetLastError();
	if ( cudaSuccess != err ){
		printf("Error! %s : %s\n",message.c_str(),cudaGetErrorString(err));
	}
}

template<typename V, typename E>
void Lprop<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);

    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    cudaMalloc(&d_adj_labels,      sizeof(int) * m);
    if(!this->myconf.onestep){
    	cudaMalloc(&d_sec_labels,         sizeof(int) * n);
    	cudaMalloc(&d_sec_adj_labels,      sizeof(int) * m);
    }
    cudaMalloc(&d_if_update,sizeof(bool) * n);
//for x_y exchange
    if(this->myconf.if_inverse){
		cudaMalloc(&d_neighborkey_out,sizeof(int) * m);
		cudaMalloc(&d_neighborval_out,sizeof(int) * m);
    }
}

template<typename V, typename E>
void Lprop<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_vertices);
    cudaFree(d_offsets);
    cudaFree(d_labels);
    cudaFree(d_sec_labels);

    cudaFree(d_counter);
    cudaFree(d_if_update);
    cudaFree(d_adj_labels);
    cudaFree(d_sec_adj_labels);
    cudaFree(d_neighborkey_out);
    cudaFree(d_neighborval_out);
}

// Return the number of labels updated
template<typename V, typename E>
int Lprop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
int Lprop<V, E>::binary_search(E left, E right,std::vector<V> offsets, std::vector<V> vertices, V target)
{
	E ans = 0;
	if(offsets[vertices[left]+1]-offsets[vertices[left]]<target){
		ans =-1;
	}else{
		if(offsets[vertices[right]+1]-offsets[vertices[right]]>target){
			ans =right;
		}else{
			while (offsets[vertices[left]+1]-offsets[vertices[left]]>target && offsets[vertices[right]+1]-offsets[vertices[right]]<target)
			{
				if(offsets[vertices[(left+right+1)/2]+1]-offsets[vertices[(left+right+1)/2]]>target){
					left = (left+right+1)/2;
					ans = left;
				}else{
					right = (left+right-1)/2;
					ans = right;
				}
			}
		}
		ans++;
	}
	if(ans != -1){
		while(offsets[vertices[ans-1]+1]-offsets[vertices[ans-1]]<=target){
			ans--;
		}
	}
	return ans;
}


template<typename V, typename E>
double Lprop<V, E>::perform_lp(V n,E m,int niter)
{
	//init flag as false
	bool flag = false;
	// init odd_even_mask as true
	bool odd_even_mask = true;
	// whether to optimate load label
	bool opt_load = this->myconf.opt_load;
	// whether one step look up, set opt_load = true before set onestep false.
	bool onestep = this->myconf.onestep;
	// whether inverse when load label.
	bool if_inverse =  this->myconf.if_inverse;
	// when to use Boolean array.
	int changeiter = this->myconf.changeiter;
	int changecount = 10;// !=0

	const int VT = 1;
	const int VT_m = 1;
	const int nthreads = 128;
    const int n_blocks = divup(n,nthreads);
    const int m_blocks = divup(m,nthreads);
    const int m_blocks_div = divup(m_blocks,VT_m);
	const int nt_big = 32*4;
	int num_huge = 0;
	int num_big = 0;
	int num_medium = 0;
	int num_small = 0;
	int num_tiny = G->n;
	// record the index in vectices array.
	int huge_index = 0;
	int big_index = 0;
	int medium_index = 0;
	int small_index = 0;

    const int layer = 4;
    const E medium = 32*layer;
    const E huge = 256*10;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);
    if(if_inverse == true){
		//d_neighborkey_out 0, 1， 2，3，...
		initialize_neighbor<<<m_blocks, nthreads>>>(d_neighborkey_out, m);
		 cudaDeviceSynchronize();
		cudaMemcpy(d_neighborval_out,d_neighbors,m*sizeof(int),cudaMemcpyDeviceToDevice);
		thrust::device_ptr<int> dev_keys_ptr(d_neighborval_out);
		thrust::device_ptr<int> dev_data_ptr(d_neighborkey_out);
//		//sort by thrust
		thrust::sort_by_key(dev_keys_ptr, dev_keys_ptr + m, dev_data_ptr);
    }

//	printf("\n");
    Timer t0;
    t0.start();
    std::vector<V> vertices;
	vertices.resize(G->n);
	 for (int i=0;i<G->n;++i) {
		vertices[i] = i;
	}

	// computer nodes size
	huge_index =  binary_search( 0,  G->n-1,G->offsets,vertices,  huge);
	printf("huge_index %d \n",huge_index);
	big_index =  binary_search( 0,  G->n-1, G->offsets,vertices,  medium);
	printf("big_index %d \n",big_index);
	medium_index = binary_search( 0,  G->n-1, G->offsets,vertices,  32);
	printf("medium_index %d \n",medium_index);
	small_index = binary_search( 0,  G->n-1, G->offsets,vertices,  16);
	printf("small_index %d \n",small_index);
//
	num_huge = huge_index+1;
	num_big = big_index-huge_index;
	num_medium = medium_index-big_index;
	num_small = small_index-medium_index;
	num_tiny = G->n-1-small_index;
	printf("num_huge %d,num_big %d,num_medium %d,num_small %d,num_tiny %d \n",num_huge,num_big,num_medium,num_small,num_tiny);
	//compute number of blocks for huge nodes
	const int nt_huge = 32*8;
	int nb_huge = 0;
	// if huge nide exists
	if(num_huge>0){
		int* nb_huge_pn;
		nb_huge_pn = &nb_huge;
		cudaMalloc(&d_block_num,      sizeof(int) * (num_huge+1));
		//compute kernel
		compute_num_blocks<nt_huge><<<divup(num_huge,nthreads),nthreads>>>(d_offsets,0,huge_index, d_block_num);
		 //exclusive scan to compute nb_huge.
		thrust::exclusive_scan(thrust::device,d_block_num,d_block_num+num_huge+1,d_block_num);
		cudaMemcpy(nb_huge_pn,d_block_num+num_huge,sizeof(int), cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		// memory assignment after compute nb_huge.
		cudaMalloc(&d_block_v,      sizeof(int) * nb_huge);
		cudaMalloc(&d_block_id,      sizeof(int) * nb_huge);
		//assign each block its vertex
		assign_blocks<<<64, 128, 0, 0>>>(d_block_num, 0,huge_index, d_block_v,d_block_id);
		// global count for across block reducing
		cudaMalloc(&d_max_count_huge, sizeof(unsigned long long)*num_huge);
		cudaMemset(d_max_count_huge, 0, sizeof(unsigned long long)*num_huge);
		printf("huge node assignment success \n");
	}
	//allocate global memory for huge and big nodes
	if(num_big > 0){
	    cudaMalloc(&(d_tables.keys),sizeof(int)*G->offsets[big_index]);
	    cudaMalloc(&(d_tables.vals),sizeof(int)*G->offsets[big_index]);
	}

	// compute tiny nodes, how many edges for each warps
	int warpnumber = 1;
	std::vector<int> warp_v;
	std::vector<int> warp_begin;
	int * d_warp_v;
	int * d_warp_begin;
	int start = small_index;
	if(small_index == -1){
		//all node are tiny
		start = 0;
	}
	warp_begin.push_back(start);
	for (int i = 0; i<= num_tiny;++i){

		if(G->offsets[small_index+ i+1]-G->offsets[start]>32){

			warpnumber++;
			for(int k =0;k<32-(G->offsets[small_index+ i]-G->offsets[start]);++k){
				warp_v.push_back(-1);
			}
			for(int k = 0; k<G->offsets[small_index+ i+1] - G->offsets[small_index+ i];++k){
				warp_v.push_back(small_index+ i);
			}
			start = small_index+ i;
			warp_begin.push_back(G->offsets[start]);
		}else{

			for(int k = 0; k<G->offsets[small_index+ i+1] - G->offsets[small_index+ i];++k){
				warp_v.push_back(small_index+ i);
			}
		}

	}
//	// mask last few threads
	while(warpnumber*32-(int)warp_v.size()>0){
		warp_v.push_back(-1);
	}

	int nt_tiny = 64;
	int nb_tiny2 = divup(warpnumber, nt_tiny/32);
	cudaMalloc (&d_warp_v,sizeof(int)*(warpnumber*32));
	cudaMalloc (&d_warp_begin,sizeof(int)*warpnumber);
	cudaMemcpy(d_warp_v,&warp_v[0],sizeof(int)*(warpnumber*32),cudaMemcpyHostToDevice);
	cudaMemcpy(d_warp_begin,&warp_begin[0],sizeof(int)*warpnumber,cudaMemcpyHostToDevice);
	printf("tiny nodes assignment success \n");


	const int nt_medium  = 32*4;
	// 1 node 1 warp
	const int nb_medium1  = divup(num_medium, nt_medium/32*VT);
	const int nt_small = 32*2;
	const int nb_small = divup(num_small, nt_small/32*VT);
	const int nb_tiny  =divup(num_tiny, nt_tiny*VT);

	std::cout<< "huge setting : " << huge<< std::endl;
	std::cout<< "medium setting : " << medium<< std::endl;
	std::cout<< "huge nodes: " << num_huge<< std::endl;
    std::cout<< "big nodes: " << num_big<< std::endl;
    std::cout<< "medium nodes: " << num_medium<< std::endl;
    std::cout<< "small nodes: " << num_small<< std::endl;
    std::cout<< "tiny nodes: " << num_tiny<< std::endl;
    std::cout<< "nb_huge: "<<nb_huge<<std::endl;
    std::cout<< "nb_big: "<<(num_big+VT-1)/VT<<std::endl;
    std::cout<< "nb_medium: "<<nb_medium1<<std::endl;
    std::cout<< "nb_small: "<<nb_small<<std::endl;
    std::cout<< "nb_tiny: "<<nb_tiny2<<std::endl;
    const int hd = G-> offsets[vertices[0]+1]- G-> offsets[vertices[0]];
    std::cout<< "highest degree: "<<hd<<std::endl;
    std::cout<< "big node edges : "<<G->offsets[big_index]-G->offsets[huge_index]<<std::endl;
	cudaMemcpy(this->d_vertices, &(vertices[0]), sizeof(V) * (G->n), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	cudaStream_t stream1, stream2, stream3;
	cudaStreamCreateWithFlags(&stream1,cudaStreamNonBlocking);
	cudaStreamCreateWithFlags(&stream2,cudaStreamNonBlocking);
	cudaStreamCreateWithFlags(&stream3,cudaStreamNonBlocking);
	float t1(0.0);
    Timer t2;
    Timer t3;
    Timer t5;
    float t4(0.0);
    float t6(0.0);
    // test edge changed number
    errorCheck("before LP!");
    cudaDeviceSynchronize();

    t0.stop();
    std::cout << " precompute: " << t0.elapsed_time()<<std::endl;
    // the first iteration
    Timer t7;
    t7.start();
    first_iteration<<<n_blocks, nthreads, 0, 0>>>(d_neighbors,d_offsets,d_labels, n);
    cudaDeviceSynchronize();
    t7.stop();
    errorCheck("start of iteration!");
    std::cout << "in interation: " << 0<< " running time: " << t7.elapsed_time() <<std::endl;
    //the other iteration
    for(int j = 1;j<niter;++j){
    	errorCheck("start of new iteration!");
    	t2.start();
    	t5.start();
    	//test
//    	auto nei = new int[G->m];
//    	cudaMemcpy(nei,d_neighbors , sizeof(E) * (G->m), cudaMemcpyDeviceToHost);
//    	for(int ni = 0; ni < 100; ni++){
//    		std::cout << "----"<<nei[ni]<<"---";
//    	}
//    	break;
    	if(flag==false){
    		gather_labels<int,int><<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
//    		gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    	}else{
    		if(onestep == true){
    			if(!if_inverse){
    				gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update);
    			}else{
    				gather_labels_inverse<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighborkey_out,d_neighborval_out, d_labels, d_adj_labels, m,d_if_update);
    			}
    		}else{

					if(odd_even_mask == true){
						// if true, generate d_adj_labels.
						if(!if_inverse){
							gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update);
						}else{

							gather_labels_inverse<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighborkey_out,d_neighborval_out, d_labels, d_adj_labels, m,d_if_update );
						}

					}else{
						//if false, generate d_sec_adj_labels.
						if(!if_inverse){

							gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_sec_labels, d_sec_adj_labels, m,d_if_update);
						}else{

							gather_labels_inverse<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighborkey_out,d_neighborval_out, d_sec_labels, d_sec_adj_labels, m,d_if_update );
						}

					}

    		}
    	}
    	cudaDeviceSynchronize();
    	t5.stop();

    	if(flag == true){
    		cudaMemset(d_if_update, 0, n*sizeof(bool));
    	}

    	errorCheck("after loading label!");
     	cudaDeviceSynchronize();
    	t3.start();
    	// reset global hash table.
		if(num_big>0){
			cudaMemset((d_tables.keys),0,sizeof(int)*G->offsets[big_index]);
			cudaMemset((d_tables.vals),0,sizeof(int)*G->offsets[big_index]);
		}
		const int buffer = 600;
		const int d = 3;
		const int w = buffer/d;
    	if(!opt_load||j<changeiter){
    		//prepare for two step look up.
			if( opt_load&&onestep== false &&j== changeiter-1){
				cudaMemcpy(d_sec_labels,d_labels,sizeof(int)*n, cudaMemcpyDeviceToDevice);
				cudaMemcpy(d_sec_adj_labels,d_adj_labels,sizeof(int)*m, cudaMemcpyDeviceToDevice);
			}

			cudaDeviceSynchronize();
			if(num_huge>0){
				 l_huge_update_syn2<32,nt_huge><<<nb_huge, nt_huge, 0,0>>>( d_neighbors,d_offsets, d_adj_labels,d_vertices,d_block_v,d_block_id,d_max_count_huge,d_tables);
				 l_sub_huge_update_syn<32,nthreads><<<divup(num_huge,nthreads),nthreads,0,0>>>(d_max_count_huge ,  d_labels ,d_counter, huge_index);
			}

			if(num_big>0){
	    	   l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(huge_index,big_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_tables);
//				l_big_update_syn3<32,VT, nt_big,buffer,w,d><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(huge_index,big_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_tables);
			}
	       if(num_medium>0){
	    	   l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
	       }
			if(num_small>0){
				l_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
			}
			if(num_tiny>0){
				l_tiny_update_syn2<<<nb_tiny2, nt_tiny,0, 0>>>(d_neighbors, d_offsets, d_adj_labels,d_labels, d_counter,d_warp_v,d_warp_begin,warpnumber);
			}

    	}else{
    		flag = true;
    		// look one step further
			if(onestep == true){
				if(num_huge>0){
					l_huge_update_syn2<32,nt_huge><<<nb_huge, nt_huge, 0,0>>>( d_neighbors,d_offsets, d_adj_labels,d_vertices,d_block_v,d_block_id,d_max_count_huge,d_tables);
					l_sub_huge_update_syn<32,nthreads><<<divup(num_huge,nthreads),nthreads,0,0>>>(d_max_count_huge ,  d_labels ,d_counter, huge_index,d_if_update);
				}
				 if(num_big>0){
					 l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(huge_index,big_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_tables,d_if_update);
				 }
				 if(num_medium>0){
					 l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_if_update);
				 }
				if(num_small>0){
					l_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_if_update);
				}
				if(num_tiny>0){
					l_tiny_update_syn2<<<nb_tiny2, nt_tiny,0, 0>>>(d_neighbors, d_offsets, d_adj_labels,d_labels, d_counter,d_warp_v,d_warp_begin,warpnumber,d_if_update);
				}

    		}else{
    			//look two steps further
    			if(odd_even_mask){
    				// first iteration using d_labels, output d_sec_labels
    				if(num_huge>0){
						l_huge_update_syn2<32,nt_huge><<<nb_huge, nt_huge, 0,0>>>( d_neighbors,d_offsets, d_adj_labels,d_vertices,d_block_v,d_block_id,d_max_count_huge,d_tables);
						l_sub_huge_update_syn<32,nthreads><<<divup(num_huge,nthreads),nthreads,0,0>>>(d_max_count_huge ,  d_sec_labels ,d_counter, huge_index,d_if_update);
    				}
    				 if(num_big>0){
    					 l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(huge_index,big_index,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter, d_tables,d_if_update);
    				 }
    				 if(num_medium>0){
    					 l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,medium,d_if_update);
    				 }
    				 if(num_small>0){
    					 l_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_if_update);
    				 }
    				 if(num_tiny>0){
    					 l_tiny_update_syn2<<<nb_tiny2, nt_tiny,0, 0>>>(d_neighbors, d_offsets, d_adj_labels,d_sec_labels, d_counter,d_warp_v,d_warp_begin,warpnumber,d_if_update);
    				 }
    			}else{
    				// second iteration using d_sec_labels, output d_labels
    				if(num_huge > 0){
    					l_huge_update_syn2<32,nt_huge><<<nb_huge, nt_huge, 0,0>>>( d_neighbors,d_offsets, d_sec_adj_labels,d_vertices,d_block_v,d_block_id,d_max_count_huge,d_tables);
    					l_sub_huge_update_syn<32,nthreads><<<divup(num_huge,nthreads),nthreads,0,0>>>(d_max_count_huge ,  d_labels ,d_counter, huge_index,d_if_update);
    				}
    				if(num_big>0){
    					l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(huge_index,big_index,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter, d_tables,d_if_update);
    				}
    				if(num_medium){
    					l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,medium,d_if_update);
    				}
    				if(num_small){
    					l_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_if_update);
    				}
    				if(num_tiny){
    					l_tiny_update_syn2<<<nb_tiny2, nt_tiny,0, 0>>>(d_neighbors, d_offsets, d_sec_adj_labels,d_labels, d_counter,d_warp_v,d_warp_begin,warpnumber,d_if_update);
    				}

    			}
    			odd_even_mask = !odd_even_mask;
    		}
		}


    	cudaDeviceSynchronize();
    	t3.stop();
    	t2.stop();
    	//printf function

        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        t6+=t5.elapsed_time();

        cudaDeviceSynchronize();
        changecount = get_count();
        std::cout << "in interation: " <<j<< ", node update : "<<changecount<< " total time(t1): "<<t2.elapsed_time()<<", processing time(t2): "<<t3.elapsed_time()<<", load label time: "<<t5.elapsed_time()<<std::endl;

        printf("\n");
    }
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    std::cout << "Mix, all time: "<< t1 + t7.elapsed_time()<< " propogation time: "<< t4+t7.elapsed_time()<< " loading label: " << t6<<std::endl;
    return t1 + t7.elapsed_time();
}

template<typename V, typename E>
double Lprop<V, E>::run(int niter)
{
	double running_time;
	this->preprocess();
	running_time = this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
	return running_time;
}

