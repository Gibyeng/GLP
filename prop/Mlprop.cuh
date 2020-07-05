// -*- coding: utf-8 -*-
// one step look up.

#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"

template<typename V, typename E>
class Mlprop {
public:
	using GraphT = CSRGraph<V, E>;
	Mlprop(GraphT* G):G(G){comparator.mp= this;}
    void run(int niter);
    int get_count();
    bool resort_compare( E i,E j);
    void testerror(int niter);
    struct Comparator {
    		Mlprop *mp;
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

    int  *d_stability; //n
    bool *d_candidate; //n
    const int mapwidth = 3200;
    int *d_gcounter; //n
    int *d_gcounter_temp; //n
    int *d_gcounterlist; //m
    int *d_temp_gcounter; //n
    bool *d_if_update_gcounter;//n
    int *d_dif_counter; //1
};

template<typename V, typename E>
void Mlprop<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);

	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);

}
template<typename V, typename E>
void Mlprop<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void Mlprop<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);
    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_temp_labels,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));

    cudaMalloc(&d_dif_counter,        sizeof(int) * 1);
    cudaMemset(d_dif_counter, 0, sizeof(int));

    cudaMalloc(&d_adj_labels,      sizeof(int) * m);

    cudaMalloc(&d_if_update,sizeof(bool) * n);
    cudaMalloc(&d_gcounter,sizeof(int) * n);
    cudaMemset(d_gcounter,0,sizeof(int) * n);
    cudaMalloc(&d_gcounter_temp,sizeof(int) * n);
    cudaMemset(d_gcounter_temp,0,sizeof(int) * n);

    cudaMalloc(&d_gcounterlist,sizeof(int) * m);
    cudaMemset(d_gcounterlist,0,sizeof(int) * m);
    cudaMalloc(&d_if_update_gcounter,sizeof(bool) * n);
    cudaMemset(d_if_update_gcounter,0,sizeof(bool) * n);


}
template<typename V, typename E>
void Mlprop<V, E>::free_gmem()
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
int Mlprop<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
void Mlprop<V, E>::perform_lp(V n,E m,int niter)
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
    float delta = 0;
    const int alpha = 0;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);
    Timer t0;
    t0.start();
    std::vector<V> vertices;
    vertices.resize(G->n);
	 for (int i=0;i<G->n;++i) {
		vertices[i] = i;
	}
	std::sort(vertices.begin(),vertices.end(),comparator);

	// computer num_medium
	E left = 0;
	E right = G->n-1;
	if(G->offsets[vertices[left]+1]-G->offsets[vertices[left]]<=32){
		num_medium =0;
	}else{
		if(G->offsets[vertices[right]+1]-G->offsets[vertices[right]]>=32){
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
//	unsigned long long* d_big_val;
//	int* d_big_key;
//	cudaMalloc(&d_big_val,      sizeof(unsigned long long) * num_big);
//	cudaMalloc(&d_big_key,      sizeof(int) * num_big);

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
    	if(flag == false){
    		// for each node update its label or global count
    		// update all its labels
    		ml_gather_labels<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
    		// update all its counts
    		ml_gather_counter<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_adj_labels, d_gcounter,d_gcounterlist, m);


    	}else{
    		// use one step look up to check difference on neigbhor labels and global counts.
    		compute_if_gcounter<<<n_blocks, nthreads>>>(d_gcounter,d_gcounter_temp,d_if_update_gcounter,alpha, n);

    		cudaDeviceSynchronize();
    		//load global label based on boolean array
    		ml_gather_labels<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m, d_if_update);
    		// load  global counter based on two boolean arrays, d_neighbors and d_if_update decise whether label changed, and d_if_update_gcounter decide whether reload.
    		ml_gather_counter2<int,int,VT_m><<<divup(m_blocks,VT_m), nthreads, 0, 0>>>(d_neighbors,d_adj_labels, d_gcounter,d_gcounterlist, m, d_if_update,d_if_update_gcounter);

    	}
    	// reset boolean array
    	if(flag == true){
    		//d_if_update_gcounter[i] = 0, not to update point i
    		cudaMemset(d_if_update_gcounter, 0, n*sizeof(bool));
    		cudaMemset(d_if_update, 0, n*sizeof(bool));
    	}
    	// copy d_gcounter to d_gcounter_temp
    	cudaMemcpy(d_gcounter_temp,d_gcounter,sizeof(int)* G->n, cudaMemcpyDeviceToDevice);
    	cudaMemset(d_gcounter, 0, n*sizeof(int));
     	cudaDeviceSynchronize();
    	t3.start();

    	if(j < 10){
    		//load counter by O(V)
//			ml_big_update_syn3<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
//			ml_medium_update_syn3<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,delta);
//			ml_small_update_syn3<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_gcounter,d_gcounter_temp,delta);
			//load counter by O(E)
    		ml_big_update_syn<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounterlist,delta);
			ml_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounterlist,delta);
			ml_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_gcounter,d_gcounterlist,delta);
    	}else{
    		flag= true;
    		//load counter by O(E) first then by O(V), compute if_update array.
    		ml_big_update_syn2<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,d_gcounterlist,delta,alpha, d_if_update);
    		ml_medium_update_syn2<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,d_gcounterlist,delta,alpha, d_if_update);
    		ml_small_update_syn2<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_gcounter,d_gcounter_temp,d_gcounterlist,delta,alpha, d_if_update);
//    		ml_big_update_syn<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounterlist,delta);
//    		ml_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounterlist,delta);
//    	    ml_small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_gcounter,d_gcounterlist,delta);
//    		ml_big_update_syn3<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(num_huge,num_huge+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,d_gcounter,d_gcounter_temp,delta);
//			ml_medium_update_syn3<32,VT,nt_medium,layer><<<nb_medium1, nt_medium,0, 0>>>(num_huge+num_big,num_huge+num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium,d_gcounter,d_gcounter_temp,delta);
//			ml_small_update_syn3<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_huge+num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_gcounter,d_gcounter_temp,delta);
    	}

    	cudaDeviceSynchronize();
//    	int *gcounter = new int [20];
//    	int *label = new int [10];
//
//    	cudaMemcpy(label,d_labels,sizeof(int)*10, cudaMemcpyDeviceToHost);
//    	cudaMemcpy(gcounter,d_gcounter_temp,sizeof(int)* 20, cudaMemcpyDeviceToHost);
//    	cudaDeviceSynchronize();
//    	std::cout << "in interation: " <<j<< std::endl;
//    	std::cout << "label: ";
//		for(int k = 0;k< 10;k++){
//			printf(" %d ",label[k]);
//
//		}
//		std::cout<<std::endl;
//    	std::cout << "gcount:";
//    	for(int k = 0;k< 20;k++){
//    		printf(" %d ",gcounter[k]);
//
//    	}
//    	std::cout<<std::endl;

    	t3.stop();
        t2.stop();
        t1+=t2.elapsed_time();
        t4+=t3.elapsed_time();
        std::cout << "in interation: " <<j<< ", update : "<<get_count()<<", total time(t1): "<<t2.elapsed_time()<<", remove random access(t2): "<<t3.elapsed_time()<<std::endl;
    }
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    std::cout << "Mix, t1: "<< t1<< " t2: "<< t4<< " t1-t2: " << t1-t4<<std::endl;
}

template<typename V, typename E>
void Mlprop<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}


// computer error ratio per niter
template<typename V, typename E>
void Mlprop<V, E>::testerror(int niter){
}
