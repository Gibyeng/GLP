// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../kernel.cuh"


template<typename V, typename E>
class HTprop {
public:
	using GraphT = CSRGraph<V, E>;
	HTprop(GraphT* G):G(G){};
    void run(int niter);
    int get_count();
    void testerror(int niter);
protected:
	void preprocess();
    void postprocess();
    void perform_lp(int n, int m,int niter, int v_offset=0);
    void init_gmem(int n, int m);
    void free_gmem();

    // Attributes
    GraphT* G;
    V* labels;
    
    int *d_neighbors;    // m
    int *d_offsets;      // n + 1
    int *d_labels;       // n
    int *d_counter;      // 1
    GlobalHT d_tables;   // m * 2
    double scale_factor;
};

template<typename V, typename E>
void HTprop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(int)* G->m, cudaMemcpyHostToDevice);
//	int a[G->m];
//	int b[G->n+1];
//	cudaMemcpy(a,d_neighbors,sizeof(int)*G->m, cudaMemcpyDeviceToHost);
//	std::cout<<"a is " <<a[0]<<std::endl;
	 //sizeof E should be 4
	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(int) * (G->n+1), cudaMemcpyHostToDevice);
//	cudaMemcpy(b,d_offsets,sizeof(int)*(G->n+1), cudaMemcpyDeviceToHost);
//	for(int j=0;j<G->n+1;j++){
//		std::cout<<"b is " <<b[j]<<" ";
//	}

}
template<typename V, typename E>
void HTprop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void HTprop<V, E>::init_gmem(int n, int m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(int) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    //scale_factor = 1.1;
    this->scale_factor = 1;
    int tm = m * scale_factor;
    // assign hash table
    cudaMalloc(&(d_tables.keys), sizeof(uint32_t) * tm);
    cudaMalloc(&(d_tables.vals), sizeof(uint32_t) * tm);

}
template<typename V, typename E>
void HTprop<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_offsets);
    cudaFree(d_labels);
    cudaFree(d_counter);

    cudaFree(d_tables.keys);
    cudaFree(d_tables.vals);
}

// Return the number of labels updated
template<typename V, typename E>
int HTprop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
void HTprop<V, E>::perform_lp(int n, int m,int niter, int v_offset)
{
    int tm = m * scale_factor;
    cudaMemset(d_tables.keys, 0, sizeof(uint32_t) * tm);
    cudaMemset(d_tables.vals, 0, sizeof(uint32_t) * tm);

    const int nthreads = 128;
    const int n_blocks = divup(n, nthreads);
    const int m_blocks = divup(m, nthreads);

    const int nt = 64;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);
    Timer t1; t1.start();
    for(int j = 0;j<niter;++j){
    	ht_update<nt, nt * 7><<<n, nt, 0>>>(d_neighbors, d_offsets, d_labels, d_tables, d_counter, 0);
    	//clear d_tables
    	cudaDeviceSynchronize();
    	cudaMemset(d_tables.keys, 0, sizeof(uint32_t) * tm);
    	cudaMemset(d_tables.vals, 0, sizeof(uint32_t) * tm);
    	cudaDeviceSynchronize();
    }
    t1.stop();
    std::cout << "hash table : "<< t1.elapsed_time()<<std::endl;
}

template<typename V, typename E>
void HTprop<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter,0);
	this->postprocess();
}

// computer error ratio per niter
template<typename V, typename E>
void HTprop<V, E>::testerror(int niter){
	this->preprocess();
	const int nthreads = 128;
	const int n_blocks = divup(G->n, nthreads);
	const int m_blocks = divup(G->m, nthreads);
	const int nt = 64;
	V* temp_labels;
	initialize_labels<<<n_blocks, nthreads>>>(d_labels, G->n);
	cudaDeviceSynchronize();
	//various of CPU classic method
	std::vector<V> vertices(G->n);
	this->labels = new V[G->n];
	auto prelabels = new V[G->n];
	 for (auto i: range(G->n)) {
		vertices[i] = i;
		labels[i] = i;
		prelabels[i] = i;
	}
	 for (auto i: range(niter)) {

	        V nupdates = 0;

	        //  copy the cpu result to gpu
	        cudaMemcpy(this->d_labels, labels,sizeof(int)* (G->n), cudaMemcpyHostToDevice);
	        cudaDeviceSynchronize();
	        // run cpu
	        for (auto i: range(G->n)) {
	        		prelabels[i] = labels[i];
	        	}
	        for (auto u: vertices) {

	            std::map<V, int> label_count;
	            V max_label = labels[u];
	            int max_count = 0;

	            for (auto v: G->iterate_neighbors(u)) {
	                //label starts with 0 but graph id starts with 1.
	                V label = labels[v];
	                int c = ++label_count[label];
	                if (max_count < c){
	                    max_count = c;
	                    max_label = label;
	                }
	            }

	            if (labels[u] != max_label) {
	                labels[u] = max_label;
	                ++nupdates;
	            }
	        }
	        // labels--cpu results, d_labels-- gpu results
	        ht_update<nt, nt * 7><<<G->n, nt, 0>>>(d_neighbors, d_offsets, d_labels, d_tables, d_counter, 0);
	        cudaDeviceSynchronize();
	        int tm = G->m * scale_factor;
			cudaMemset(d_tables.keys, 0, sizeof(uint32_t) * tm);
			cudaMemset(d_tables.vals, 0, sizeof(uint32_t) * tm);
			cudaDeviceSynchronize();
	        //copy results back to CPU
	        temp_labels = new V [G->n];

	        cudaMemcpy(temp_labels,d_labels,sizeof(int)* G->n, cudaMemcpyDeviceToHost);

	        cudaDeviceSynchronize();
	        int dif_count = 0;
	        for(int j = 0;j<G->n;j++){
	        	if(labels[j]!=temp_labels[j]){
	        		++dif_count;
	        	}
	        }
            std::cout << "in iteration "<< i<< ", label error is "<< dif_count<< " error ratio is "<< (double) dif_count/G->n<< std::endl;
//	        std::cout << "GPU ";
//	        std::cout << "in iteration "<< i<< " different is " ;
//	        for(int i = 0;i<5;i++){
//	        	std::cout<<temp_labels[i]<<" ";
//	        }
//	        std::cout<<std::endl;
//	        std::cout << "CPU ";
//			std::cout << "in iteration "<< i<< " different is " ;
//			for(int i = 0;i<5;i++){
//				std::cout<<labels[i]<<" ";
//			}
//			std::cout<<std::endl;
	        delete[] temp_labels;
	    }
	 delete[] labels;
	 delete[] prelabels;
	this->postprocess();
}
