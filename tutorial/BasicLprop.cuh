// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
template<typename V, typename E>
class BasicLprop {
public:
	using GraphT = CSRGraph<V, E>;
	BasicLprop(GraphT* G):G(G){}
    void run(int niter);
    int get_count();

    void testerror(int niter);
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
    int *d_counter;      // 1
    double scale_factor;
    int *d_adj_labels; //m
    bool *d_if_update; //n

    // for huge nodes

    int *d_block_num;
	int* d_block_v;
	int* d_block_id;
	GlobalHT d_tables;

	//for glolbal reduce
	int  *d_max_count_huge;
    // for y-x sort
    int *d_neighborkey_out; //m
    int *d_neighborval_out; //m
};

template<typename V, typename E>
void BasicLprop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);

	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);


}
template<typename V, typename E>
void BasicLprop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void BasicLprop<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);
    cudaMalloc(&d_adj_labels,         sizeof(int) *m);
    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    cudaMalloc(&d_max_count_huge,sizeof(int)*n);
    cudaMalloc(&d_tables.keys,sizeof(int)*m);
    cudaMalloc(&d_tables.vals,sizeof(int)*m);

    cudaMalloc(&d_neighborkey_out,sizeof(int) * m);
    cudaMalloc(&d_neighborval_out,sizeof(int) * m);
}
template<typename V, typename E>
void BasicLprop<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_vertices);
    cudaFree(d_offsets);
    cudaFree(d_labels);

    cudaFree(d_counter);
    cudaFree(d_adj_labels);

    cudaFree(d_neighborkey_out);
    cudaFree(d_neighborval_out);
}

// Return the number of labels updated
template<typename V, typename E>
int BasicLprop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
int BasicLprop<V, E>::binary_search(E left, E right,std::vector<V> offsets, std::vector<V> vertices, V target)
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
void BasicLprop<V, E>::perform_lp(V n,E m,int niter)
{

	const int VT = 7;
	const int VT_m = 1;
	const int nthreads = 128;
    const int n_blocks = divup(n,nthreads);
    const int m_blocks = divup(m,nthreads);
    const int m_blocks_div = divup(m_blocks,VT_m);



    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);

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

	//compute number of blocks for nodes
	const int nt_huge = 32*2;
	int nb_huge = 0;
	int* nb_huge_pn;
	nb_huge_pn = &nb_huge;
	cudaMalloc(&d_block_num,      sizeof(int) * (G->n+1));
	//compute kernel
	compute_num_blocks<nt_huge><<<divup(G->n,nthreads),nthreads>>>(d_offsets,0,G->n-1, d_block_num);
	 //exclusive scan to compute nb_huge.
	thrust::exclusive_scan(thrust::device,d_block_num,d_block_num+G->n+1,d_block_num);
	cudaMemcpy(nb_huge_pn,d_block_num+G->n,sizeof(int), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	// memory assignment after compute nb_huge.
	cudaMalloc(&d_block_v,      sizeof(int) * nb_huge);
	cudaMalloc(&d_block_id,      sizeof(int) * nb_huge);
	//assign each block its vertex
	assign_blocks<<<64, 128, 0, 0>>>(d_block_num, 0,G->n-1, d_block_v,d_block_id);

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
    std::cout << " precompute: " << t0.elapsed_time()<<std::endl;

    for(int j = 0;j<niter;++j){

    	t2.start();
    	t5.start();


    	gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);

    	cudaDeviceSynchronize();
    	t5.stop();

     	cudaDeviceSynchronize();
    	t3.start();
    	cudaMemset(d_tables.keys,0,sizeof(int)*G->m);
    	cudaMemset(d_tables.vals,0,sizeof(int)*G->m);
    	cudaMemset(d_max_count_huge,0,sizeof(int)*G->n);
    	cudaDeviceSynchronize();

//		l_huge_update_syn2<32,nt_huge><<<nb_huge, nt_huge, 0,0>>>(d_neighbors,d_offsets, d_adj_labels,d_labels,d_vertices,d_block_v,d_block_id,d_max_count_huge, d_counter,d_tables);
//		l_sub_huge_update_syn<32,nthreads><<<n_blocks,nthreads,0,0>>>(d_max_count_huge ,  d_labels ,d_counter, n);
    	cudaDeviceSynchronize();
    	t3.stop();
    	t2.stop();
    	//printf function
		int *label1 = new int [8];
//		int *label2 = new int [20];
//		int *gcount = new int [20];
		cudaMemcpy(label1,d_labels,sizeof(int)*8, cudaMemcpyDeviceToHost);
//		cudaMemcpy(label2,d_sec_labels,sizeof(int)*20, cudaMemcpyDeviceToHost);
//		cudaMemcpy(gcount,d_gcounter,sizeof(int)*20, cudaMemcpyDeviceToHost);
//		cudaDeviceSynchronize();
		std::cout << "in interation: " <<j<< std::endl;
		std::cout << "label: ";
		for(int k = 0;k< 8;k++){
			printf(" %d ",label1[k]);
		}
		std::cout<<std::endl;

        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        t6+=t5.elapsed_time();

        cudaDeviceSynchronize();
        std::cout << "in interation: " <<j<< ", node update : "<<get_count()<< " total time(t1): "<<t2.elapsed_time()<<", processing time(t2): "<<t3.elapsed_time()<<", load label time: "<<t5.elapsed_time()<<std::endl;
//
//        printf("\n");
    }
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    std::cout << "Mix, all time: "<< t1<< " propogation time: "<< t4<< " loading label: " << t6<<std::endl;
}

template<typename V, typename E>
void BasicLprop<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}


// computer error ratio per niter
template<typename V, typename E>
void BasicLprop<V, E>::testerror(int niter){

}
