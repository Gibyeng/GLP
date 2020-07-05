// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
template<typename V, typename E>
class Ml2prop {
public:
	using GraphT = CSRGraph<V, E>;
	Ml2prop(GraphT* G):G(G){comparator.mp= this;}
    void run(int niter);
    int get_count();
    bool resort_compare( E i,E j);
    void testerror(int niter);
    int binary_search(E left,E right,std::vector<V>  offsets,std::vector<V> vertices,V target);
    struct Comparator {
    		Ml2prop *mp;
    		bool operator () (E i,E j) {
    			//std::cout << "i : "<< i << " j :"<< j<< " offset i:"<<mp->G->offsets[i+1]-mp->G->offsets[i]<< " offset j:"<<mp->G->offsets[j+1]-mp->G->offsets[j]<< std::endl;
    			return mp->G->offsets[j+1]-mp->G->offsets[j] < mp->G->offsets[i+1]-mp->G->offsets[i];
    		}
    	}comparator;

    GraphT* G;
protected:
    void preprocess();
    void postprocess();
    void perform_lp(V n,E m,int niter);
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
    double scale_factor;
    int *d_adj_labels; //m
    int *d_sec_adj_labels;//m
    bool *d_if_update; //n
    int *d_gcounter; //n
    int *d_gcounter_temp; //n

//    GlobalHT d_tables;

    int *d_neighborkey_out; //m
    int *d_neighborval_out; //m
};

template<typename V, typename E>
void Ml2prop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);

	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);


}
template<typename V, typename E>
void Ml2prop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void Ml2prop<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);
    cudaMalloc(&d_sec_labels,         sizeof(int) * n);
    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    cudaMalloc(&d_adj_labels,      sizeof(int) * m);
    cudaMalloc(&d_sec_adj_labels,      sizeof(int) * m);
    cudaMalloc(&d_if_update,sizeof(bool) * n);

    cudaMalloc(&d_gcounter,sizeof(int) * n);
    cudaMalloc(&d_gcounter_temp,sizeof(int) * n);
    cudaMemset(d_gcounter,0,sizeof(int) * n);
    cudaMemset(d_gcounter_temp,0,sizeof(int) * n);
    cudaMalloc(&d_neighborkey_out,sizeof(int) * m);
    cudaMalloc(&d_neighborval_out,sizeof(int) * m);
}
template<typename V, typename E>
void Ml2prop<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_vertices);
    cudaFree(d_offsets);
    cudaFree(d_labels);
    cudaFree(d_sec_labels);
    cudaFree(d_gcounter);
    cudaFree(d_gcounter_temp);
    cudaFree(d_counter);
    cudaFree(d_if_update);
    cudaFree(d_adj_labels);
    cudaFree(d_sec_adj_labels);
    cudaFree(d_neighborkey_out);
    cudaFree(d_neighborval_out);
}

// Return the number of labels updated
template<typename V, typename E>
int Ml2prop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
int Ml2prop<V, E>::binary_search(E left, E right,std::vector<V> offsets, std::vector<V> vertices, V target)
{
	E ans = 0;
	if(offsets[vertices[left]+1]-offsets[vertices[left]]<target){
		ans =0;
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
	while(offsets[vertices[ans-1]+1]-offsets[vertices[ans-1]]<=target){
		ans--;
	}
	return ans;
}


template<typename V, typename E>
void Ml2prop<V, E>::perform_lp(V n,E m,int niter)
{
	//init flag as false
	bool flag = false;
	// init odd_even_mask as true
	bool odd_even_mask = true;
	// whether to optimate load label
	bool opt_load = false;
	// whether one step look up, set opt_load = true before set onestep false.
	bool onestep = true;
	// whether inverse when load label.
	bool if_inverse = false;
	// when to use Boolean array.
	int changeiter = 10;
	// delta is for layered LP
	float delta = 0;
	const int VT = 7;
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

    const int k =1;
    const int layer = 4;
    const E medium = 32*layer;
    const E huge = 100000;
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
    //end sort
//    cudaDeviceSynchronize();
//	int *neighborval_out = new int [20];
//	int *neighborkey_out = new int [20];
//	cudaMemcpy(neighborval_out,d_neighborval_out,sizeof(int)*20, cudaMemcpyDeviceToHost);
//	cudaMemcpy(neighborkey_out,d_neighborkey_out,sizeof(int)*20, cudaMemcpyDeviceToHost);
//	 cudaDeviceSynchronize();
//	 std::cout<< std::endl<< "val: ";
//	for(int k = 0;k< 20;k++){
//		printf(" %d ",neighborval_out[k]);
//	}
//	std::cout<< std::endl<< "key: ";
//	for(int k = 0;k< 20;k++){
//		printf(" %d ",neighborkey_out[k]);
//	}

//	printf("\n");
    Timer t0;
    t0.start();
    std::vector<V> vertices;
	vertices.resize(G->n);
	 for (int i=0;i<G->n;++i) {
		vertices[i] = i;
	}
	std::sort(vertices.begin(),vertices.end(),comparator);

	// computer
//	huge_index =  binary_search( 0,  G->n-1,G->offsets,vertices,  huge);
	printf("huge_index %d \n",huge_index);
	big_index =  binary_search( 0,  G->n-1, G->offsets,vertices,  medium);
	printf("big_index %d \n",big_index);
	medium_index = binary_search( 0,  G->n-1, G->offsets,vertices,  32);
	printf("medium_index %d \n",medium_index);
	small_index = binary_search( 0,  G->n-1, G->offsets,vertices,  16);
	printf("small_index %d \n",small_index);

	num_huge = huge_index+1;
	num_big = big_index-huge_index;
	num_medium = medium_index-big_index;
	num_small = small_index-medium_index;
	num_tiny = G->n-small_index;
	const int nt_medium  = 32*2;
	// 1 node 1 warp
	const int nb_medium1  = divup(num_medium, nt_medium/32*VT);
	// 1 node 1 block
	const int nb_medium2  = divup(num_medium, nt_medium);
	const int nt_small = 32*2;
	const int nb_small = divup(num_small, nt_small/32*VT);
	const int nt_tiny  = 32*2;
	const int nb_tiny  =divup(num_tiny, nt_tiny*VT);
	const int nt_tiny2  = 32*2;
	const int nb_tiny2  =divup(G->n, nt_tiny*VT);
	std::cout<< "huge setting : " << huge<< std::endl;
	std::cout<< "medium setting : " << medium<< std::endl;
	std::cout<< "huge nodes: " << num_huge<< std::endl;
    std::cout<< "big nodes: " << num_big<< std::endl;
    std::cout<< "medium nodes: " << num_medium<< std::endl;
    std::cout<< "small nodes: " << num_small<< std::endl;
    std::cout<< "tiny nodes: " << num_tiny<< std::endl;
    const int hd = G-> offsets[vertices[0]+1]- G-> offsets[vertices[0]];
    std::cout<< "highest degree: "<<hd<<std::endl;
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

    cudaDeviceSynchronize();

    t0.stop();
    std::cout << "sorted and precompute: " << t0.elapsed_time()<<std::endl;

    for(int j = 0;j<niter;++j){

    	t2.start();
    	t5.start();

    	if(flag==false){
//    		gather_labels<int,int><<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    		gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
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
    	cudaMemcpy(d_gcounter_temp,d_gcounter,sizeof(int)* G->n, cudaMemcpyDeviceToDevice);
    	cudaMemset(d_gcounter, 0, n*sizeof(int));
    	if(flag == true){
    		cudaMemset(d_if_update, 0, n*sizeof(bool));

    	}


     	cudaDeviceSynchronize();
    	t3.start();

    	if(!opt_load||j<changeiter){
    		//prepare for two step look up.
			if( opt_load&&onestep== false &&j== changeiter-1){
				cudaMemcpy(d_sec_labels,d_labels,sizeof(int)*n, cudaMemcpyDeviceToDevice);
				cudaMemcpy(d_sec_adj_labels,d_adj_labels,sizeof(int)*m, cudaMemcpyDeviceToDevice);
			}
			ml_big_update_syn3<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,big_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
			ml_medium_update_syn3<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,delta);
			ml_small_update_syn3<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
//			ml_tiny_update_syn3_cmp<VT><<<nb_tiny, nt_tiny,0, 0>>>(nb_small,G->n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
			ml_tiny_update_syn3<VT><<<nb_tiny, nt_tiny,0, 0>>>(small_index,G->n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
			ml_tiny_update_syn3_cmp<VT><<<nb_tiny, nt_tiny,0, 0>>>(small_index,G->n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
//			ml_tiny_update_syn4<VT><<<nb_tiny2, nt_tiny2,0, 0>>>(0,G->n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
    	}else{
    		flag = true;
    		// look one step further
			if(onestep == true){

				ml_big_update_syn3<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,big_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp, delta,d_if_update);
				ml_medium_update_syn3<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,delta,d_if_update);
				ml_small_update_syn3<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
//				ml_tiny_update_syn3<VT><<<nb_tiny, nt_tiny,0, 0>>>(small_index,G->n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
				ml_tiny_update_syn4<VT><<<nb_tiny2, nt_tiny2,0, 0>>>(0,G->n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    		}else{
    			//look two steps further
    			if(odd_even_mask){
    				// first iteration using d_labels, output d_sec_labels
    				ml_big_update_syn3<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,big_index,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    				ml_medium_update_syn3<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,delta,d_if_update);
    				ml_small_update_syn3<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
//    				ml_tiny_update_syn3<VT><<<nb_tiny, nt_tiny,0, 0>>>(small_index,G->n,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    				ml_tiny_update_syn4<VT><<<nb_tiny2, nt_tiny2,0, 0>>>(0,G->n,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    			}else{
    				// second iteration using d_sec_labels, output d_labels
    				ml_big_update_syn3<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,big_index,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    				ml_medium_update_syn3<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,delta,d_if_update);
    				ml_small_update_syn3<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,small_index,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
//    				ml_tiny_update_syn3<VT><<<nb_tiny, nt_tiny,0, 0>>>(small_index,G->n,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    				ml_tiny_update_syn4<VT><<<nb_tiny2, nt_tiny2,0, 0>>>(0,G->n,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta,d_if_update);
    			}
    			odd_even_mask = !odd_even_mask;
    		}
		}
    	cudaDeviceSynchronize();
    	t3.stop();
    	t2.stop();
    	//printf function
//		int *label1 = new int [20];
//		int *label2 = new int [20];
//		int *gcount = new int [20];
//		cudaMemcpy(label1,d_labels,sizeof(int)*20, cudaMemcpyDeviceToHost);
//		cudaMemcpy(label2,d_sec_labels,sizeof(int)*20, cudaMemcpyDeviceToHost);
//		cudaMemcpy(gcount,d_gcounter,sizeof(int)*20, cudaMemcpyDeviceToHost);
//		cudaDeviceSynchronize();
//		std::cout << "in interation: " <<j<< std::endl;
//		std::cout << "label: ";
//		for(int k = 0;k< 20;k++){
//			printf(" %d ",label1[k]);
//		}
//		std::cout<<std::endl;
//		std::cout << "d_label_sec:";
//		for(int k = 0;k< 20;k++){
//			printf(" %d ",label2[k]);
//
//		}
//		std::cout<<std::endl;
//		std::cout << "gcount: ";
//		for(int k = 0;k< 20;k++){
//			printf(" %d ",gcount[k]);
//		}
//		std::cout<<std::endl;
//
//    	cudaDeviceSynchronize();
//
        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        t6+=t5.elapsed_time();

        cudaDeviceSynchronize();
        std::cout << "in interation: " <<j<< ", node update : "<<get_count()<< " total time(t1): "<<t2.elapsed_time()<<", processing time(t2): "<<t3.elapsed_time()<<", load label time: "<<t5.elapsed_time()<<std::endl;

        printf("\n");
    }
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    std::cout << "Mix, all time: "<< t1<< " propogation time: "<< t4<< " loading label: " << t6<<std::endl;
}

template<typename V, typename E>
void Ml2prop<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}


// computer error ratio per niter
template<typename V, typename E>
void Ml2prop<V, E>::testerror(int niter){

}
