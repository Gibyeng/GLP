// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
template<typename V, typename E>
class Base4prop {
public:
	using GraphT = CSRGraph<V, E>;
	Base4prop(GraphT* G):G(G){}
    double run(int niter);
    int get_count();

    void experiment(const int niter,const bool opt_load, const bool onestep, const bool if_inverse, std::string filename);
    int binary_search(E left,E right,std::vector<V>  offsets,std::vector<V> vertices,V target);

    GraphT* G;
protected:
    void preprocess();
    void postprocess();
    double perform_lp(V n,E m,int niter);
    void init_gmem(V n,E m);
    void free_gmem();
    void errorCheck(std::string message);
    // Attributes
    V* labels;
    int *d_vertices;     // n
    int *d_neighbors;    // m
    E *d_offsets;      // n + 1
    int *d_labels;       // n

    int *d_counter;      // 1

    int *d_adj_labels; //m

    bool *d_if_update; //n

    // for huge nodes
    int *d_block_num;
	int* d_block_v;
	int* d_block_id;
	//for medium nodes
	int *d_warp_num;
	int* d_warp_v;
	int* d_warp_id;
	GlobalHT d_tables; //m*2
	//for glolbal reduce
	unsigned long long  *d_max_count_huge;

    //
    int* d_gcounter_r; // n
	int* d_gcounter_w; // n
};

template<typename V, typename E>
void Base4prop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);
}
template<typename V, typename E>
void Base4prop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void Base4prop<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);

    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    cudaMalloc(&d_adj_labels,      sizeof(int) * m);

    cudaMalloc(&d_if_update,sizeof(bool) * n);

;
    cudaMalloc(&d_gcounter_r,sizeof(int) * n);
    cudaMalloc(&d_gcounter_w,sizeof(int) * n);
    cudaMemset(d_gcounter_r,0,sizeof(int)*n);
    cudaMemset(d_gcounter_w,0,sizeof(int)*n);
}
template<typename V, typename E>
void Base4prop<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_vertices);
    cudaFree(d_offsets);
    cudaFree(d_labels);


    cudaFree(d_counter);
    cudaFree(d_if_update);
    cudaFree(d_adj_labels);


}

template<class V, class E>
void Base4prop<V, E>::errorCheck(std::string message){
	auto err = cudaGetLastError();
	if ( cudaSuccess != err ){
		printf("Error! %s : %s\n",message.c_str(),cudaGetErrorString(err));
	}
}

// Return the number of labels updated
template<typename V, typename E>
int Base4prop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
int Base4prop<V, E>::binary_search(E left, E right,std::vector<V> offsets, std::vector<V> vertices, V target)
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
double Base4prop<V, E>::perform_lp(V n,E m,int niter)
{
	const int VT = 7;
	const int VT_m = 1;
	const int nthreads = 128;
    const int n_blocks = divup(n,nthreads);
    const int m_blocks = divup(m,nthreads);
    const int m_blocks_div = divup(m_blocks,VT_m);
	const int nt_big = 32*4;

	int d_iter = 10;
	int delta = 1;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);


//	printf("\n");
    Timer t0;
    t0.start();
    std::vector<V> vertices;
	vertices.resize(G->n);
	 for (int i=0;i<G->n;++i) {
		vertices[i] = i;
	}

	//allocate global memory for nodes

	cudaMalloc(&(d_tables.keys),sizeof(int)*G->m);
	cudaMalloc(&(d_tables.vals),sizeof(int)*G->m);



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


    t0.stop();
    std::cout << " precompute: " << t0.elapsed_time()<<std::endl;
    // the first iteration
    Timer t7;
    t7.start();
//    first_iteration<<<n_blocks, nthreads, 0, 0>>>(d_neighbors,d_offsets,d_labels, n);
    cudaDeviceSynchronize();
    t7.stop();
    errorCheck("before LP!");
    for(int s = 0; s < d_iter; s++){
    	delta = 2*delta;
    	int delta_star = delta +1;
    	if(s!=0){
    		initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);
    	}
		//the other iteration
		for(int j = 0;j<niter;++j){

			t2.start();
			t5.start();
			gather_labels<int,int,VT_m><<<m_blocks_div, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
			 errorCheck("before gahter label!");
			cudaDeviceSynchronize();
			t5.stop();
			t3.start();
			// reset global hash table.
			cudaMemset((d_tables.keys),0,sizeof(int)*m);
			cudaMemset((d_tables.vals),0,sizeof(int)*m);
			cudaDeviceSynchronize();
			// update global counter;
			if(j!=0){
				cudaMemcpy(d_gcounter_r, d_gcounter_w,sizeof(int)*n,cudaMemcpyDeviceToDevice);
				cudaMemset(d_gcounter_w, 0, sizeof(int)*n);
			}
			ml_big_update_syn_no_shared<32,VT, nt_big><<<(n+VT-1)/VT, nt_big, 0,0>>>(0,m-1,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_tables,d_gcounter_r,d_gcounter_w, delta_star);
			errorCheck("before update!");
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
    }
    std::cout << "Mix, all time: "<< t1 + t7.elapsed_time()<< " propogation time: "<< t4+t7.elapsed_time()<< " loading label: " << t6<<std::endl;
    return t1 + t7.elapsed_time();
}

template<typename V, typename E>
double Base4prop<V, E>::run(int niter)
{
	double running_time =0;
	this->preprocess();
	running_time = this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
	return running_time;
}

