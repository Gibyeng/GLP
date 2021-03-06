// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
template<typename V, typename E>
class MBase2prop {
public:
	using GraphT = CSRGraph<V, E>;
	MBase2prop(GraphT* G):G(G){}
    void run(int niter);
    int get_count();

    void experiment(const int niter,const bool opt_load, const bool onestep, const bool if_inverse, std::string filename);
    int binary_search(E left,E right,std::vector<V>  offsets,std::vector<V> vertices,V target);

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
};

template<typename V, typename E>
void MBase2prop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);
}
template<typename V, typename E>
void MBase2prop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void MBase2prop<V, E>::init_gmem(V n,E m)
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
void MBase2prop<V, E>::free_gmem()
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
int MBase2prop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
int MBase2prop<V, E>::binary_search(E left, E right,std::vector<V> offsets, std::vector<V> vertices, V target)
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
void MBase2prop<V, E>::perform_lp(V n,E m,int niter)
{
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
	int num_small = G->n;
	// record the index in vectices array.
	int huge_index = 0;
	int big_index = 0;
	int medium_index = 0;
	int small_index = 0;

    const int layer = 4;
    const E medium = 32*layer;
    const E huge = 10000;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);

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

//
	num_huge = huge_index+1;
	num_big = big_index-huge_index;
	num_medium = medium_index-big_index;
	num_small = G->n-1-medium_index;
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
	}
	//allocate global memory for huge and big nodes
	if(num_big > 0){
	    cudaMalloc(&(d_tables.keys),sizeof(int)*G->offsets[big_index]);
	    cudaMalloc(&(d_tables.vals),sizeof(int)*G->offsets[big_index]);
	}
	const int nt_medium  = 32*4;
	// 1 node 1 warp
	const int nb_medium1  = divup(num_medium, nt_medium/32*VT);
	const int nt_small = 32*2;
	const int nb_small = divup(num_small, nt_small/32*VT);


	std::cout<< "huge setting : " << huge<< std::endl;
	std::cout<< "medium setting : " << medium<< std::endl;
	std::cout<< "huge nodes: " << num_huge<< std::endl;
    std::cout<< "big nodes: " << num_big<< std::endl;
    std::cout<< "medium nodes: " << num_medium<< std::endl;
    std::cout<< "nb_big: "<<(num_big+VT-1)/VT<<std::endl;
    std::cout<< "nb_medium: "<<nb_medium1<<std::endl;
    std::cout<< "nb_small: "<<nb_small<<std::endl;
    const int hd = G-> offsets[vertices[0]+1]- G-> offsets[vertices[0]];
    std::cout<< "highest degree: "<<hd<<std::endl;
	cudaMemcpy(this->d_vertices, &(vertices[0]), sizeof(V) * (G->n), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	float t1(0.0);
    Timer t2;
    Timer t3;
    Timer t5;
    float t4(0.0);
    float t6(0.0);
    // test edge changed number

    cudaDeviceSynchronize();

    t0.stop();
    std::cout << " precompute: " << t0.elapsed_time()<<std::endl;
    // the first iteration
    Timer t7;
    t7.start();
    first_iteration<<<n_blocks, nthreads, 0, 0>>>(d_neighbors,d_offsets,d_labels, n);
    cudaDeviceSynchronize();
    t7.stop();
    std::cout << "in interation: " << 0<< " running time: " << t7.elapsed_time() <<std::endl;
    //the other iteration
    for(int j = 1;j<niter;++j){
    	t2.start();
    	t5.start();
		gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    	cudaDeviceSynchronize();
    	t5.stop();
    	t3.start();
    	// reset global hash table.
		if(num_big>0){
			cudaMemset((d_tables.keys),0,sizeof(int)*G->offsets[big_index]);
			cudaMemset((d_tables.vals),0,sizeof(int)*G->offsets[big_index]);
		}
		if(num_huge>0){
			 l_huge_update_syn2<32,nt_huge><<<nb_huge, nt_huge, 0,0>>>( d_neighbors,d_offsets, d_adj_labels,d_vertices,d_block_v,d_block_id,d_max_count_huge,d_tables);
			 l_sub_huge_update_syn<32,nthreads><<<divup(num_huge,nthreads),nthreads,0,0>>>(d_max_count_huge ,  d_labels ,d_counter, huge_index);
		}
		l_big_update_syn2<32,VT, nt_big,600><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(huge_index,big_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_tables);
		l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(big_index,medium_index,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
		l_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(medium_index,m,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);

    	cudaDeviceSynchronize();
    	t3.stop();
    	t2.stop();
    	//printf function

        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        t6+=t5.elapsed_time();

        cudaDeviceSynchronize();
        int changecount = get_count();
        std::cout << "in interation: " <<j<< ", node update : "<<changecount<< " total time(t1): "<<t2.elapsed_time()<<", processing time(t2): "<<t3.elapsed_time()<<", load label time: "<<t5.elapsed_time()<<std::endl;

        printf("\n");
    }

    std::cout << "Mix, all time: "<< t1 + t7.elapsed_time()<< " propogation time: "<< t4+t7.elapsed_time()<< " loading label: " << t6<<std::endl;
}

template<typename V, typename E>
void MBase2prop<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}

