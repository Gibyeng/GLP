// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"

template<typename V, typename E>
class Mixprop_p {
public:
	using GraphT = CSRGraph<V, E>;
	Mixprop_p(GraphT* G):G(G){comparator.mp= this;}
    void run(int niter);
    int get_count();
    bool resort_compare( E i,E j);
    void testerror(int niter);
    struct Comparator {
    		Mixprop_p *mp;
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
    int *d_temp_labels;       // n
    int *d_counter;      // 1
    double *d_prob;
    double *d_rand_prob;
    double scale_factor;
    bool *d_if_update;
    int *d_adj_labels;
};

template<typename V, typename E>
void Mixprop_p<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);

}
template<typename V, typename E>
void Mixprop_p<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void Mixprop_p<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);
    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_temp_labels,         sizeof(int) * n);
    cudaMalloc(&d_prob,         sizeof(double) * n);
    cudaMalloc(&d_rand_prob,         sizeof(double) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));
    cudaMemset(d_prob,0,sizeof(double));
    cudaMalloc(&d_if_update,sizeof(bool) * n);
    cudaMalloc(&d_adj_labels,      sizeof(int) * m);

}
template<typename V, typename E>
void Mixprop_p<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_vertices);
    cudaFree(d_offsets);
    cudaFree(d_labels);
    cudaFree(d_temp_labels);
    cudaFree(d_counter);
    cudaFree(d_if_update);
    cudaFree(d_adj_labels);

}

// Return the number of labels updated
template<typename V, typename E>
int Mixprop_p<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
void Mixprop_p<V, E>::perform_lp(V n,E m,int niter)
{
	const int VT = 5;
	const int nthreads = 128;
    const int n_blocks = divup(n,nthreads);
	const int m_blocks = divup(G->m, nthreads);
	const int nt_big = 32*4;
	int num_big = n*1;
	int num_medium = n*0;
	int num_small = G->n - num_big - num_medium;
    const int k =1;
    const int d = 5;
    const E medium = 32*d;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);
    std::vector<V> vertices;
    vertices.resize(G->n);
    for (int i=0;i<G->n;++i) {
   		vertices[i] = i;
   	}
	Timer t6;
	t6.start();
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
	}
	num_small = G->n-num_medium;
	num_medium = num_medium-num_big;
	t6.stop();
	std::cout << "resort cost: "<<t6.elapsed_time()<<std::endl;
	const int nt_medium  = 32*2;
	const int nb_medium1  = divup(num_medium, nt_medium /32);
	const int nb_medium2  = divup(num_medium, nt_medium);
	const int nt_small = 32*2;
	const int nb_small = divup(num_small, nt_small/32);
	std::cout<< "medium setting : " << medium<< std::endl;
    std::cout<< "big nodes: " << num_big<< std::endl;
    std::cout<< "medium nodes: " << num_medium<< std::endl;
    std::cout<< "small nodes: " << num_small<< std::endl;
	cudaMemcpy(this->d_vertices, &(vertices[0]), sizeof(int) * (G->n), cudaMemcpyHostToDevice);
	cudaStream_t stream1, stream2, stream3;
	cudaStreamCreate (&stream1);
	cudaStreamCreate (&stream2);
	cudaStreamCreate (&stream3);
	float t1(0.0);
    Timer t2;
    Timer t3;
    Timer t7;
    float t4(0.0);
    t7.start();
	gen_prob<<<n_blocks, nthreads>>>(d_rand_prob, n);
	cudaDeviceSynchronize();
	t7.stop();
	std::cout << "random cost: "<<t7.elapsed_time()<<std::endl;

    for(int j = 0;j<niter;++j){
    	t2.start();
    	if(j < 4){
    		gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    	}else{
//    		gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update);
    	}
    	cudaMemset(d_if_update, 0, n*sizeof(bool));
		cudaDeviceSynchronize();
		t3.start();
//    	big_update_syn<32, nt_big,32,200,k><<<(num_big+k-1)/k, nt_big, 0,stream1>>>(0,num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter,d_prob,d_rand_prob);
//    	//big_update_syn<32, nt_big,32,200,k><<<(num_big+k-1)/k, nt_big, 0,stream1>>>(0,num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
//    	medium_update_syn<32,nt_medium><<<nb_medium1, nt_medium,0, stream2>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter,d_prob,d_rand_prob);
//    	small_update_syn<32,nt_small><<<nb_small, nt_small,0, stream3>>>(num_big+num_medium,n,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter,d_prob,d_rand_prob);
    	if(j < 3){
    		big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_prob,d_rand_prob);
    		//medium_update_syn4<32,nt_medium,d><<<nb_medium1, nt_medium,0, 0>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_prob,d_rand_prob);
    		small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_prob,d_rand_prob);
    	}
    	else{
			big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update,d_prob,d_rand_prob);
			//medium_update_syn4<32,nt_medium,d><<<nb_medium1, nt_medium,0, 0>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium, d_if_update,d_prob,d_rand_prob);
			small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update,d_prob,d_rand_prob);
		}
    	cudaDeviceSynchronize();
		t3.stop();
		t2.stop();
		t1+=t2.elapsed_time();
		t4+=t3.elapsed_time();
		std::cout << "in interation: " <<j<< ", update : "<<get_count()<<", total time(t1): "<<t2.elapsed_time()<<", remove random access(t2): "<<t3.elapsed_time()<<std::endl;
	}

	std::cout << "Mix, t1: "<< t1<< " t2: "<< t4<<std::endl;
}

template<typename V, typename E>
void Mixprop_p<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}


