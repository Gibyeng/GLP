// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"

template<typename V, typename E>
class Mixprop_g {
public:
	using GraphT = CSRGraph<V, E>;
	Mixprop_g(GraphT* G):G(G){comparator.mp= this;}
    void run(int niter);
    int get_count();
    bool resort_compare( E i,E j);
    void testerror(int niter);
    struct Comparator {
    	    Mixprop_g *mp;
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
    SharedCM<200,3>* d_cm;
    double scale_factor;
    int *d_adj_labels; //m
};

template<typename V, typename E>
void Mixprop_g<V,E>::preprocess(){
	this->init_gmem(G->n, G->m);
	cudaMemcpy(this->d_neighbors, &(G->neighbors[0]),sizeof(V)*G->m, cudaMemcpyHostToDevice);
    //int a[G->m];
	//int64_t b[G->n+1];
//	int* a = new int [G->m];
//	int64_t* b =new int64_t [G->n+1];
//	cudaMemcpy(a,d_neighbors,sizeof(V)*G->m, cudaMemcpyDeviceToHost);
//	for(int j=0;j<30;j++){
//		std::cout<<"a is " <<a[j] <<std::endl;
//	}
	cudaMemcpy(this->d_offsets, &(G->offsets[0]), sizeof(E) * (G->n+1), cudaMemcpyHostToDevice);
//	cudaMemcpy(b,d_offsets,sizeof(E)*(G->n+1), cudaMemcpyDeviceToHost);
//	for(int j=0;j<10;j++){
//		std::cout<<"b is " <<b[j]<<" ";
//	}
//	std::cout<<"size of " <<G->offsets.size()<<" ";
}
template<typename V, typename E>
void Mixprop_g<V, E>::postprocess()
{
    this->free_gmem();
}

template<typename V, typename E>
void Mixprop_g<V, E>::init_gmem(V n,E m)
{
    cudaMalloc(&d_neighbors,      sizeof(int) * m);
    cudaMalloc(&d_offsets,        sizeof(E) * (n + 1));
    cudaMalloc(&d_labels,         sizeof(int) * n);
    cudaMalloc(&d_vertices,         sizeof(int) * n);
    cudaMalloc(&d_temp_labels,         sizeof(int) * n);
    cudaMalloc(&d_counter,        sizeof(int) * 1);
    cudaMemset(d_counter, 0, sizeof(int));
    cudaMalloc(&d_cm, sizeof(SharedCM<200,3>)* n);
    cudaMemset(d_cm, 0, sizeof(SharedCM<200,3>)* n);
    cudaMalloc(&d_adj_labels,      sizeof(int) * m);
}
template<typename V, typename E>
void Mixprop_g<V, E>::free_gmem()
{
    cudaFree(d_neighbors);
    cudaFree(d_vertices);
    cudaFree(d_offsets);
    cudaFree(d_labels);
    cudaFree(d_temp_labels);
    cudaFree(d_counter);


}

// Return the number of labels updated
template<typename V, typename E>
int Mixprop_g<V, E>::get_count()
{
    int counter;
    cudaMemcpy(&counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter, 0, sizeof(int));

    return counter;
}

template<typename V, typename E>
void Mixprop_g<V, E>::perform_lp(V n,E m,int niter)
{
	const int VT = 5;
	const int nthreads = 128;
	const int n_blocks = divup(n,nthreads);
	const int m_blocks = divup(m,nthreads*1	);
	const int nt_big = 32*4;
	int num_big = n*1;
	int num_medium = n*0;
	int num_small = G->n - num_big - num_medium;
	const int k =1;
	const int d = 10;
	const E medium = 32*d;
    initialize_labels<<<n_blocks, nthreads>>>(d_labels, n);
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
	const int nt_medium  = 32*2;
	const int nb_medium1  = divup(num_medium, nt_medium /32);
	const int nb_medium2  = divup(num_medium, nt_medium);
	const int nt_small = 32*2;
	const int nb_small = divup(num_small, nt_small/32);
	std::cout<< "medium setting : " << medium<< std::endl;
    std::cout<< "big nodes: " << num_big<< std::endl;
    std::cout<< "medium nodes: " << num_medium<< std::endl;
    std::cout<< "small nodes: " << num_small<< std::endl;
	cudaMemcpy(this->d_vertices, &(vertices[0]), sizeof(V) * (G->n), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	cudaStream_t stream1, stream2, stream3;
	//	cudaStreamCreate (&stream1);
	//	cudaStreamCreate (&stream2);
	//	cudaStreamCreate (&stream3);
	cudaStreamCreateWithFlags(&stream1,cudaStreamNonBlocking);
	cudaStreamCreateWithFlags(&stream2,cudaStreamNonBlocking);
	cudaStreamCreateWithFlags(&stream3,cudaStreamNonBlocking);
    Timer t1;
    Timer t2;
    Timer t3;
    float t4(0.0);
    t1.start();
    for(int j = 0;j<niter;++j){
    	t2.start();
		gather_labels<<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
		cudaDeviceSynchronize();
		t3.start();
    	big_update_syn3_g<32,VT, nt_big,32,256><<<(num_big+VT-1)/VT, nt_big, 0,0>>>(0,num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_cm );
    	//big_update_syn3_g<32, nt_big,32,200,k><<<(num_big+k-1)/k, nt_big, 0,stream1>>>(0,num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter, d_cm );
    	//big_update_syn<32, nt_big,32,200,k><<<(num_big+k-1)/k, nt_big, 0,stream1>>>(0,num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter );
    	//big_update_syn3<32, nt_big,32,200><<<num_big, nt_big, 0,stream1>>>(0,num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
    	//medium_update_syn<32,nt_medium><<<nb_medium1, nt_medium,0, stream2>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
    	//medium_update_syn2<32,nt_medium><<<nb_medium2, nt_medium,0, stream2>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
    	//medium_update_syn3<32,nt_medium><<<nb_medium2, nt_medium,0, stream2>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_labels,d_labels,d_vertices, d_counter);
    	//medium_update_syn5<32,nt_medium,d><<<nb_medium1, nt_medium,0, 0>>>(num_big,num_medium+num_big,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter,medium);
    	small_update_syn<32,VT,nt_small><<<nb_small, nt_small,0, 0>>>(num_big+num_medium,n,d_neighbors, d_offsets, d_adj_labels,d_labels,d_vertices, d_counter);
    	cudaDeviceSynchronize();
    	t3.stop();
    	t2.stop();
    	t4+=t3.elapsed_time();
    	std::cout << "in interation: " <<j<< ", update : "<<get_count()<<", total time(t1): "<<t2.elapsed_time()<<", remove random access(t2): "<<t3.elapsed_time()<<std::endl;
    }
    t1.stop();
    //std::cout << "in interation: " <<niter<< " updated : "<<this->get_count()<<std::endl;
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    std::cout << "global Mix, t1: "<< t1.elapsed_time()<< " t2: "<< t4<<std::endl;
}

template<typename V, typename E>
void Mixprop_g<V, E>::run(int niter)
{
	this->preprocess();
	this->perform_lp( G->n,  G->m, niter);
	this->postprocess();
}


// computer error ratio per niter
template<typename V, typename E>
void Mixprop_g<V, E>::testerror(int niter){
	this->preprocess();
	const int nthreads = 128;
	const int n_blocks = divup(G->n, nthreads);
	const int m_blocks = divup(G->m, nthreads);
    const int nt = 64;
    const int n_b = divup(G->n, nt/32);
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
	 std::cout << "before sort "<< vertices[0]<<" ,"<< vertices[1]<<std::endl;
	 std::sort(vertices.begin(),vertices.end(),comparator);
	 std::cout << "after sort "<< vertices[0]<<" ,"<< vertices[1]<<std::endl;
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
	        // labels--cpu results, d_labels-- gpu results
	        big_update_syn_g<32, nt,32,200,1><<<G->n, nt, 0,0>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter,d_cm);
	        //big_update_syn<32, nt,32,200,1><<<G->n, nt, 0,0>>>(0,G->n,d_neighbors, d_offsets, d_labels,d_temp_labels,d_vertices, d_counter);
	        cudaDeviceSynchronize();
	        //copy d_temp_labels to d_labels
			cudaMemcpy(d_labels,d_temp_labels,sizeof(int)* G->n, cudaMemcpyDeviceToDevice);
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
	        std::cout << "GPU ";
	        std::cout << "in iteration "<< i<< " different is " ;
	        for(int i = 0;i<36;i++){
	        	std::cout<<temp_labels[i]<<" ";
	        }
	        std::cout<<std::endl;
	        std::cout << "CPU ";
			std::cout << "in iteration "<< i<< " different is " ;
			for(int i = 0;i<36;i++){
				std::cout<<labels[i]<<" ";
			}
			std::cout<<std::endl;
	        delete[] temp_labels;
	    }
	 delete[] labels;
	 delete[] prelabels;
	this->postprocess();
}
