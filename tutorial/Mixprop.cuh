// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

template<typename V, typename E>
class Mixprop {
public:
	using GraphT = CSRGraph<V, E>;
	Mixprop(GraphT* G):G(G){comparator.mp= this;}
    void run(int niter);
    int get_count();
    bool resort_compare( E i,E j);
    void testerror(int niter);
    struct Comparator {
    	    Mixprop *mp;
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

//    GlobalHT d_tables;

    int *d_neighborkey_out; //m
    int *d_neighborval_out; //m
};

template<typename V, typename E>
void Mixprop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);

	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);


}
template<typename V, typename E>
void Mixprop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void Mixprop<V, E>::init_gmem(V n,E m)
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

    cudaMalloc(&d_neighborkey_out,sizeof(int) * m);
    cudaMalloc(&d_neighborval_out,sizeof(int) * m);
}
template<typename V, typename E>
void Mixprop<V, E>::free_gmem()
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
int Mixprop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
void Mixprop<V, E>::perform_lp(V n,E m,int niter)
{
	//init flag as false
	bool flag = false;
	// init odd_even_mask as true
	bool odd_even_mask = true;
	// whether one step look up
	bool onestep = false;
	// whether inverse when load label.
	bool if_inverse = true;
	// when to use Boolean array.
	int changeiter = 10;

	const int VT = 7;
	const int VT_m = 1;
	const int nthreads = 128;
    const int n_blocks = divup(n,nthreads);
    const int m_blocks = divup(m,nthreads);
	const int nt_big = 32*4;
	int num_huge = 0;
	int num_big = 0;
	int num_medium = 0;
	int num_small = G->n - num_huge - num_big - num_medium;
    const int k =1;
    const int layer = 5;
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
   		//sort by thrust
   		thrust::sort_by_key(dev_keys_ptr, dev_keys_ptr + m, dev_data_ptr);

   }
    //end sort

	printf("\n");
    Timer t0;
    t0.start();
    std::vector<V> vertices;
    vertices.resize(G->n);
	 for (int i=0;i<G->n;++i) {
		vertices[i] = i;
	}
	std::sort(vertices.begin(),vertices.end(),comparator);

	// computer num_big
	E left = 0;
	E right = G->n-1;
	if(G->offsets[vertices[left]+1]-G->offsets[vertices[left]]<32){
		num_medium =0;
	}else{
		if(G->offsets[vertices[right]+1]-G->offsets[vertices[right]]>32){
			num_medium =n;
		}else{
			while (G->offsets[vertices[left]+1]-G->offsets[vertices[left]]>32 && G->offsets[vertices[right]+1]-G->offsets[vertices[right]]<32)
			{
				if(G->offsets[vertices[(left+right+1)/2]+1]-G->offsets[vertices[(left+right+1)/2]]>32){
					left = (left+right+1)/2;
					num_medium = left;
				}else{
					right = (left+right-1)/2;
					num_medium = right;
				}
			}
		}
		num_medium++;
	}
	while(G->offsets[vertices[num_medium-1]+1]-G->offsets[vertices[num_medium-1]]<=32){
		num_medium--;
	}
	left = 0;
	right = G->n-1;
	if(G->offsets[vertices[left]+1]-G->offsets[vertices[left]]<=medium){
		num_big =0;
	}else{
		if(G->offsets[vertices[right]+1]-G->offsets[vertices[right]]>=medium){
			num_big =n;
		}else{
			while (G->offsets[vertices[left]+1]-G->offsets[vertices[left]]>medium && G->offsets[vertices[right]+1]-G->offsets[vertices[right]]<medium)
			{
				if(G->offsets[vertices[(left+right+1)/2]+1]-G->offsets[vertices[(left+right+1)/2]]>medium){
					left = (left+right+1)/2;
					num_big = left;

				}else{
					right = (left+right-1)/2;
					num_big = right;
				}
			}
		}
		num_big++;
	}
	while(G->offsets[vertices[num_big-1]+1]-G->offsets[vertices[num_big-1]]<=medium){
		num_big--;
	}
	// computer num_huge
	 left = 0;
	 right = G->n-1;
	if(G->offsets[vertices[left]+1]-G->offsets[vertices[left]]<huge){
		num_huge =0;
	}else{
		if(G->offsets[vertices[right]+1]-G->offsets[vertices[right]]>huge){
			num_huge =n;
		}else{
			while (G->offsets[vertices[left]+1]-G->offsets[vertices[left]]>huge && G->offsets[vertices[right]+1]-G->offsets[vertices[right]]<huge)
			{
				if(G->offsets[vertices[(left+right+1)/2]+1]-G->offsets[vertices[(left+right+1)/2]]>huge){
					left = (left+right+1)/2;
					num_huge = left;
				}else{
					right = (left+right-1)/2;
					num_huge = right;
				}
			}
		}
		num_huge++;
	}
	num_small = G->n-num_medium;
	num_medium = num_medium-num_big;
	num_big = num_big - num_huge;
	const int nt_medium  = 32*2;
	// 1 node 1 warp
	const int nb_medium1  = divup(num_medium, nt_medium/32*VT);
	// 1 node 1 block
	const int nb_medium2  = divup(num_medium, nt_medium);
	const int nt_small = 32*2;
	const int nb_small = divup(num_small, nt_small/32*VT);
	std::cout<< "huge setting : " << huge<< std::endl;
	std::cout<< "medium setting : " << medium<< std::endl;
	std::cout<< "huge nodes: " << num_huge<< std::endl;
    std::cout<< "big nodes: " << num_big<< std::endl;
    std::cout<< "medium nodes: " << num_medium<< std::endl;
    std::cout<< "small nodes: " << num_small<< std::endl;
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
    float t4(0.0);


    cudaDeviceSynchronize();

    t0.stop();
    std::cout << "sorted and precompute: " << t0.elapsed_time()<<std::endl;

    for(int j = 0;j<niter;++j){
    	t2.start();


    	if(flag==false){
//    		gather_labels<int,int><<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    		gather_labels<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);

    	}else{
    		if(onestep == true){
    			if(!if_inverse){
    				gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update);
    			}else{
    				gather_labels_inverse<<<m_blocks, nthreads, 0, 0>>>(d_neighborkey_out,d_neighborval_out, d_labels, d_adj_labels, m,d_if_update);
    			}
    		}else{
    			if(odd_even_mask == true){
    				// if true, generate d_adj_labels.
    				if(!if_inverse){
    					gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update);
    				}else{
    					gather_labels_inverse<<<m_blocks, nthreads, 0, 0>>>(d_neighborkey_out,d_neighborval_out, d_labels, d_adj_labels, m,d_if_update );
    				}

				}else{
					//if false, generate d_sec_adj_labels.
					if(!if_inverse){
						gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_sec_labels, d_sec_adj_labels, m,d_if_update);
					}else{
						gather_labels_inverse<<<m_blocks, nthreads, 0, 0>>>(d_neighborkey_out,d_neighborval_out, d_sec_labels, d_sec_adj_labels, m,d_if_update );
					}

    			}
    		}

    	}

    	if(flag == true){
    		cudaMemset(d_if_update, 0, n*sizeof(bool));
    	}


     	cudaDeviceSynchronize();
    	t3.start();

    	if(j<changeiter){
    		//prepare for two step look up.
			if(onestep== false &&j== changeiter-1){
				cudaMemcpy(d_sec_labels,d_labels,sizeof(int)*n, cudaMemcpyDeviceToDevice);
				cudaMemcpy(d_sec_adj_labels,d_adj_labels,sizeof(int)*m, cudaMemcpyDeviceToDevice);
			}
			big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
			medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);

			small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);

    	}else{
    		flag = true;
    		// look one step further
			if(onestep == true){
				big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update);
				medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium, d_if_update);
				small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update);
    		}else{
    			//look two steps further
    			if(odd_even_mask){
    				// first iteration using d_labels, output d_sec_labels
    				big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_if_update);
					medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,medium,d_if_update);
					small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_sec_labels,d_vertices, d_counter,d_if_update);
    			}else{
    				// second iteration using d_sec_labels, output d_labels
					big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_if_update);
					medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,medium,d_if_update);
					small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_sec_adj_labels,d_labels,d_vertices, d_counter,d_if_update);
    			}
    			odd_even_mask = !odd_even_mask;
    		}
		}
    	cudaDeviceSynchronize();
		int *label1 = new int [20];
		int *label2 = new int [20];
		cudaMemcpy(label1,d_labels,sizeof(int)*20, cudaMemcpyDeviceToHost);
		cudaMemcpy(label2,d_sec_labels,sizeof(int)*20, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		std::cout << "in interation: " <<j<< std::endl;
		std::cout << "label: ";
		for(int k = 0;k< 20;k++){
			printf(" %d ",label1[k]);
		}
		std::cout<<std::endl;
		std::cout << "d_label_sec:";
		for(int k = 0;k< 20;k++){
			printf(" %d ",label2[k]);

		}
		std::cout<<std::endl;

    	cudaDeviceSynchronize();
    	t3.stop();
        t2.stop();
        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        std::cout << "in interation: " <<j<< ", update : "<<get_count()<<", total time(t1): "<<t2.elapsed_time()<<", remove random access(t2): "<<t3.elapsed_time()<<std::endl;

        printf("\n");
    }
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    std::cout << "Mix, t1: "<< t1<< " t2: "<< t4<< " t1-t2: " << t1-t4<<std::endl;
}

template<typename V, typename E>
void Mixprop<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}


// computer error ratio per niter
template<typename V, typename E>
void Mixprop<V, E>::testerror(int niter){

}
