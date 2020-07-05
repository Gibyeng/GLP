// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <cub/cub.cuh>

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
    int *d_temp_labels;       // n
    int *d_counter;      // 1
    double scale_factor;
    int *d_adj_labels; //m
    bool *d_if_update; //n
    bool *d_if_update2; //n
    bool *d_if_update_bitmap;//n*ratio
    int  *d_stability; //n
    bool *d_candidate; //n
    const int mapwidth = 3200;
//    GlobalHT d_tables;
    int *d_neighborkey_in; //m
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
    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_temp_labels,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    cudaMalloc(&d_adj_labels,      sizeof(int) * m);

    cudaMalloc(&d_if_update,sizeof(bool) * n);
    cudaMalloc(&d_if_update_bitmap,sizeof(bool) * mapwidth);
//    cudaMalloc(&d_if_update2,sizeof(bool) * n);
//    cudaMalloc(&d_stability,sizeof(int) * n);
//    cudaMalloc(&d_candidate,sizeof(bool) * n);
//    cudaMalloc(&(d_tables.keys), sizeof(uint32_t) * m);
//    cudaMalloc(&(d_tables.vals), sizeof(uint32_t) * m);

    cudaMalloc(&d_neighborkey_in,sizeof(int) * m);
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
    cudaFree(d_temp_labels);
    cudaFree(d_counter);
    cudaFree(d_if_update);
    cudaFree(d_adj_labels);
    cudaFree(d_neighborkey_in);
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
    initialize_neighbor<<<m_blocks, nthreads>>>(d_neighborkey_in, m);
    // sort by cub
    void     *d_temp_storage = NULL;
    size_t   temp_storage_bytes = 0;
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
    		d_neighbors, d_neighborval_out,d_neighborkey_in, d_neighborkey_out, m);
    // Allocate temporary storage
    cudaMalloc(&d_temp_storage, temp_storage_bytes);
    // Run sorting operation
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
    		d_neighbors, d_neighborval_out, d_neighborkey_in, d_neighborkey_out, m);
    //end sort
    // check sorted array
    int* c_labels = new int[20];
	cudaMemcpy(c_labels,d_neighborval_out,20*sizeof(int),cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	for (int j =0;j< 20;j++){
		printf(" %d ",c_labels[j]);
	}
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

    const int mapwidth = this->mapwidth;

    cudaDeviceSynchronize();

    t0.stop();
    std::cout << "sorted and precompute: " << t0.elapsed_time()<<std::endl;
    bool flag = false;
    for(int j = 0;j<niter;++j){
    	t2.start();
    	if(flag==false){
//    		gather_labels<int,int><<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    		gather_labels<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);

    	}else{
//    		gather_labels<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    		gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update);
//    		gather_labels<int,int,VT_m,128,3200><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m,d_if_update_bitmap);
    	}
//    	cudaMemset(d_candidate, 0, n*sizeof(bool));
//    	cudaDeviceSynchronize();
//    	compute_candidate<V,E,nthreads><<<G->n, nthreads, 0, 0>>>(d_neighbors,d_offsets,d_if_update,d_candidate);
//    	cudaDeviceSynchronize();
//    	if(j >= 4){
//    		cudaMemcpy(d_if_update2, d_if_update,sizeof(bool)* (G->n), cudaMemcpyDeviceToDevice);
//    	}
    	if(flag == true){
    		cudaMemset(d_if_update, 0, n*sizeof(bool));
 //   		cudaMemset(d_if_update_bitmap, 0, mapwidth*sizeof(bool));
    	}


     	cudaDeviceSynchronize();
    	t3.start();

    	if(j<5){

//    		big_update_syn5<32, nt_big,32,256,3><<<num_big*3, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_big_val,d_big_key);
//   		//big_update_syn<32,VT, nt_big,32,200><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
//			big_update_syn2<32, nt_big,32,200><<<num_big, nt_big, 0,stream1>>>(0,num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
			big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
//medium_update_syn<32,nt_medium><<<nb_medium1, nt_medium,0, stream2>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
//			//medium_update_syn5<32,nt_medium,d><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_if_update);
//medium_update_syn2<32,nt_medium><<<nb_medium2, nt_medium,0, stream2>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
//medium_update_syn3<32,nt_medium><<<nb_medium2, nt_medium,0, stream2>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
//  		   medium_update_syn10<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
			medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
//   		medium_update_syn8<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
//    		medium_update_syn9<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
//			medium_update_syn6<32,nt_medium,256><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
			small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
    	}else{
    		flag = true;
			big_update_syn4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update);
//			big_update_syn_bitmap4<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update_bitmap,mapwidth);
//			medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium, d_if_update);
			medium_update_syn4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium, d_if_update);
//			medium_update_syn_bitmap4<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium, d_if_update_bitmap,mapwidth);
//			medium_update_syn8<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium, d_if_update);
//			medium_update_syn6<32,nt_medium,256><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update);
			small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update);
//			small_update_syn_bitmap<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_if_update_bitmap,mapwidth);
		}
    	//base line
//        cudaMemset(d_tables.keys, 0, sizeof(uint32_t) * m);
//        cudaMemset(d_tables.vals, 0, sizeof(uint32_t) * m);
//    	cudaDeviceSynchronize();
//    	update_lockfree<64><<<n, 64, 0, 0>>>( d_neighbors, d_offsets,d_adj_labels, d_labels, d_tables, d_vertices,d_counter);
//        update_lockfree<64><<<num_big, 64, 0, 0>>>(0, num_big, d_neighbors, d_offsets,d_adj_labels, d_labels, d_tables, d_vertices,d_counter);
//        update_lockfree<64><<<num_medium, 64, 0, 0>>>(num_big, num_medium+num_big, d_neighbors, d_offsets,d_adj_labels, d_labels, d_tables, d_vertices,d_counter);
//        update_lockfree<64><<<num_small, 64, 0, 0>>>(num_big+num_medium,n, d_neighbors, d_offsets,d_adj_labels, d_labels, d_tables, d_vertices,d_counter);
//    	big_update_syn5<32, nt_big,32,256,3><<<num_big*3, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_big_val,d_big_key);
//    	cudaDeviceSynchronize();
    	//muti-block
//    	copy_arr_to_label<<<divup(num_big,128),128,0,0>>>(d_big_key,d_big_val,d_labels, num_big );
//    	cudaDeviceSynchronize();
//    	for(int i = 0;i<3;i++){
//    		big_update_syn6<32, nt_big,32,256,3><<<num_big, nt_big, 0,0>>>(0,num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,i);
//    	}
    	cudaDeviceSynchronize();
    	t3.stop();
        t2.stop();
        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        std::cout << "in interation: " <<j<< ", update : "<<get_count()<<", total time(t1): "<<t2.elapsed_time()<<", remove random access(t2): "<<t3.elapsed_time()<<std::endl;
//        std::cout << "in interation: "<< j<<std::endl;
//        int* c_labels = new int[n];
//        cudaMemcpy(c_labels,d_labels,n*sizeof(int),cudaMemcpyDeviceToHost);
//        cudaDeviceSynchronize();
//        for (int j =0;j< 20;j++){
//			printf(" %d ",c_labels[j]);
//		}
//        delete[] c_labels;
//        int* c_labels = new int[12];
//        cudaMemcpy(c_labels,&d_adj_labels[G->offsets[0]],12*sizeof(int),cudaMemcpyDeviceToHost);
//        cudaDeviceSynchronize();
//        for (int j =0;j< 12;j++){
//			printf(" %d ",c_labels[j]);
//		}
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
	this->preprocess();
	const int nthreads = 128;
	const int n_blocks = divup(G->n, nthreads);
	const int m_blocks = divup(G->m, nthreads);
	const int VT = 7;
    const int nt = 64;
    const int n_b = divup(G->n, nt);
    const int n_b2 = divup(G->n, nt/32);
	const int n_b3 = divup(G->n, nt/32*VT);
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
	 //std::sort(vertices.begin(),vertices.end(),comparator);

	 cudaMemcpy(this->d_vertices, &(vertices[0]), sizeof(int) * (G->n), cudaMemcpyHostToDevice);
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
	            V max_label = prelabels[u];
	            int max_count = 0;

	            for (auto v: G->iterate_neighbors(u)) {
	                //label starts with 0 but graph id starts with 1.

	                V label = prelabels[v];
	                int c = ++label_count[label];
	                if (max_count < c){
	                    max_count = c;
	                    max_label = label;
	                }
	            }

	            if (prelabels[u] != max_label) {
	                labels[u] = max_label;
	                ++nupdates;
	            }
	        }


			cudaMemset(d_if_update, 0, G->n*sizeof(bool));
			cudaDeviceSynchronize();
	        // labels--cpu results, d_labels-- gpu results
	       //big_update_syn<32, nt,32,200,1><<<G->n, nt>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
	        //big_update_syn3<32,1, nt,32,200><<<G->n, nt>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
	        //medium_update_syn<32, nt><<<n_b2, nt>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
	        //medium_update_syn2<32, nt><<<n_b, nt>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
	        //medium_update_syn3<32, nt><<<n_b, nt>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
	        //small_update_syn<32, nt><<<n_b2, nt>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
//	        medium_update_syn4<32,nt,2><<<n_b2, nt,0, 0>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter,64);
			//medium_update_syn8_test<32,VT,nt,2><<<n_b3, nt,0, 0>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter,64);
			//update_lockfree<64><<<G->n, 64, 0, 0>>>(d_neighbors, d_offsets,d_labels, d_temp_labels, d_tables, d_vertices,d_counter);


			cudaDeviceSynchronize();
	        auto err = cudaGetLastError();
	        if(err!=cudaSuccess){printf("Error: %s \n",cudaGetErrorString(err));}
	        cudaDeviceSynchronize();
	        //copy d_temp_labels to d_labels
			cudaMemcpy(d_labels,d_temp_labels,sizeof(int)* G->n, cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
	        //copy results back to CPU
	        V* temp_labels = new V [G->n];

	    	bool* if_update = new bool [G->n];
	        cudaMemcpy(temp_labels,d_labels,sizeof(int)* G->n, cudaMemcpyDeviceToHost);
	        cudaMemcpy(if_update,d_if_update,sizeof(bool)* G->n, cudaMemcpyDeviceToHost);


	        cudaDeviceSynchronize();
	        int dif_count = 0;
	        for(int j = 0;j<G->n;j++){
	        	if(labels[j]!=temp_labels[j]){
	        		++dif_count;
	        	}
	        }
	        std::cout << "in iteration "<< i<< ", label error is "<< dif_count<< " error ratio is "<< (double) dif_count/G->n<< std::endl;

	        std::cout << "GPU ";
	        std::cout << "in iteration "<< i<< " label list: " ;
	        for(int i = 0;i<146;i++){
	        	std::cout<<temp_labels[i]<<" ";
	        }
	        std::cout<<std::endl;
	        std::cout << "CPU ";
			std::cout << "in iteration "<< i<< " label list:  " ;
			for(int i = 0;i<146;i++){
				std::cout<<labels[i]<<" ";
			}

			std::cout<<std::endl;
			std::cout<< "updated: "<< get_count()<<std::endl;
	        delete[] temp_labels;
	    }
	 delete[] labels;
	 delete[] prelabels;
	this->postprocess();
}
