// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <thrust/sort.h>
#include <vector>
#include <thrust/execution_policy.h>
template<typename V, typename E>
class MultiGprop {
public:
	using GraphT = CSRGraph<V, E>;
	MultiGprop(GraphT* G):G(G){}
    double run(int niter);
    int get_count();
    void errorCheck(std::string message);
    int binary_search(E left,E right,std::vector<V>  offsets,std::vector<V> vertices,V target);
    void GPUquery();
    GraphT* G;
protected:
    void preprocess();
    void postprocess();
    double perform_lp(V n,E m,int niter);
    void init_gmem(V n,E m);
    void free_gmem();
    // Attributes
    int n_A;
    int n_B;
    int m_A;
    int m_B;
	std::vector<int> neighborA;
	std::vector<int> neighborB;
	std::vector<int> offsetA;
	std::vector<int> offsetB;
    V* labels_A;
    V* labels_B;
    int *d_vertices_A;     // n
    int *d_vertices_B;
    int *d_neighbors_A;    // m
    int *d_neighbors_B;    // m
    E *d_offsets_A;      // n + 1
    E *d_offsets_B;
    int *d_labels_A;       // n
    int *d_labels_B;
    int* d_labels_odd_A;
    int* d_labels_even_A;
    int* d_labels_odd_B;
    int* d_labels_even_B;
    int *d_adj_labels_A; //m
    int *d_adj_labels_B; //m
    int* d_sec_labels_even_A; //n
    int* d_sec_labels_odd_A;
    int* d_sec_labels_even_B;
    int* d_sec_labels_odd_B;
    int* d_sec_adj_labels_A;//m
    int* d_sec_adj_labels_B;//m
    bool *d_if_update_even_A; //n
    bool *d_if_update_even_B; //n
    bool *d_if_update_odd_A; //n
    bool *d_if_update_odd_B;
    bool *d_if_update_B;
    int *d_counter_A;
    int *d_counter_B;
    int* d_r_offsets_A;
    int* d_r_offsets_B;
    // for huge nodes
    int *d_block_num_A;
    int *d_block_num_B;
	int* d_block_v_A;
	int* d_block_id_A;
	int* d_block_v_B;
	int* d_block_id_B;

	//for medium nodes
	int *d_warp_num_A;
	int* d_warp_v_A;
	int* d_warp_id_A;
	int *d_warp_num_B;
	int* d_warp_v_B;
	int* d_warp_id_B;
	GlobalHT d_tables_A; //big and huge edges *2
	GlobalHT d_tables_B;
	//for glolbal reduce
	unsigned long long  *d_max_count_huge_A;
	unsigned long long  *d_max_count_huge_B;
    // for y-x sort
    int *d_neighborkey_out_A; //m
    int *d_neighborkey_out_B; //m
    int *d_neighborval_out_A; //m
    int *d_neighborval_out_B; //m

    struct configure{
    	bool opt_load = false;
		// whether one step look up, set opt_load = true before set onestep false.
		bool onestep = true;
		// whether inverse when load label.
		bool if_inverse = false;
		// when to use Boolean array.
		int changeiter = 8;
    } myconf;
};

template<typename V, typename E>
int MultiGprop<V, E>::get_count()
{
    int counter_A;
    int counter_B;
    cudaMemcpy(&counter_A, d_counter_A, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&counter_B, d_counter_B, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemset(d_counter_A, 0, sizeof(int));
    cudaMemset(d_counter_B, 0, sizeof(int));
    return counter_A+counter_B;
}

template<typename V, typename E>
int MultiGprop<V, E>::binary_search(E left, E right,std::vector<V> offsets, std::vector<V> vertices, V target)
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
template<class V, class E>
void MultiGprop<V, E>::errorCheck(std::string message){
	auto err = cudaGetLastError();
	if ( cudaSuccess != err ){
		printf("Error! %s : %s\n",message.c_str(),cudaGetErrorString(err));
	}
}


template<typename V, typename E>
void MultiGprop<V, E>::GPUquery(){
	int ngpus;
	cudaGetDeviceCount(&ngpus);
	for (int i =0;i<ngpus;i++){

		cudaDeviceProp devProp;
		cudaGetDeviceProperties(&devProp, i);
		printf("Device %d has compute capability %d.%d.\n", i, devProp.major,
		devProp.minor);
	}
	int canAccessPeer;
	cudaDeviceCanAccessPeer(&canAccessPeer,1,0);
	if(canAccessPeer){
		printf("peer access is allowed\n");
	}
	cudaSetDevice(0);
	cudaDeviceEnablePeerAccess(1,0);
	cudaSetDevice(1);
	cudaDeviceEnablePeerAccess(0,0);
	errorCheck("EnablePeerAccess");
}

template<typename V, typename E>
void MultiGprop<V,E>::preprocess(){
	//compute segments for each GPU
	//odd node to GPU1, even nodes to GPU0
	int n = G->n;
	int m = G->m;
	int topA = 0;
	int topB = 0;
	offsetA.push_back(topA);
	offsetB.push_back(topB);
	for(int i =0;i<n;i++){
		int begin = G->offsets[i];
		int end = G->offsets[i+1];
		int length = end -begin;
		if(i%2==0){
			topA += length;
			offsetA.push_back(topA);
			for(int j = 0;j<length;j++){
				neighborA.push_back(G->neighbors[begin+j]);
			}
		}else{
			topB += length;
			offsetB.push_back(topB);
			for(int j = 0;j<length;j++){
				neighborB.push_back(G->neighbors[begin+j]);
			}
		}
	}
	this->n_A = offsetA.size()-1;
	this->n_B = offsetB.size()-1;
	this->m_A = neighborA.size();
	this->m_B = neighborB.size();
	printf("G : G->offsets[2] %d,G->offsets[1] %d, G->offsets[0] %d \n",G->offsets[2],G->offsets[1], G->offsets[0]);
	printf("A : offsetA[2] %d,offsetA[1] %d, offsetA[0] %d \n",offsetA[2],offsetA[1], offsetA[0]);

	// assign Gmem;
	cudaSetDevice(0);
	cudaMalloc(&d_neighbors_A,      sizeof(int) * m_A);
	cudaMalloc(&d_offsets_A,        sizeof(E) * (n_A + 1));

	cudaMalloc(&d_labels_A,         sizeof(int) * n);
	cudaMalloc(&d_labels_odd_A,         sizeof(int) * n_B);
	cudaMalloc(&d_labels_even_A,         sizeof(int) * n_A);
	cudaMalloc(&d_sec_labels_even_A,         sizeof(int) * n_A);
	cudaMalloc(&d_sec_labels_odd_A,         sizeof(int) * n_B);
	cudaMalloc(&d_vertices_A,         sizeof(int) * n);
	cudaMalloc(&d_counter_A,        sizeof(int) * 1);
	cudaMemset(d_counter_A, 0, sizeof(int));

	cudaMalloc(&d_adj_labels_A,      sizeof(int) * m_A);
	cudaMalloc(&d_sec_adj_labels_A,      sizeof(int) * m_A);
	cudaMalloc(&d_if_update_even_A,sizeof(bool) * n_A);
	cudaMalloc(&d_if_update_odd_A,sizeof(bool) * n_B);

	cudaMalloc(&d_r_offsets_A,        sizeof(E) * (n_A+1));
	if(this->myconf.if_inverse==true){
		cudaMalloc(&d_neighborkey_out_A,sizeof(int) * m_A);
		cudaMalloc(&d_neighborval_out_A,sizeof(int) * m_A);
	}

	cudaSetDevice(1);
	cudaMalloc(&d_neighbors_B,      sizeof(int) * (m_B));
	cudaMalloc(&d_offsets_B,        sizeof(E) * (n_B + 1));

	cudaMalloc(&d_labels_B,         sizeof(int) * n);
	cudaMalloc(&d_labels_odd_B,         sizeof(int) * n_B);
	cudaMalloc(&d_labels_even_B,         sizeof(int) * n_A);
	cudaMalloc(&d_sec_labels_even_B,        sizeof(int) * n_A);
	cudaMalloc(&d_sec_labels_odd_B,         sizeof(int) * n_B);
	cudaMalloc(&d_vertices_B,         sizeof(int) * n);
	cudaMalloc(&d_counter_B,        sizeof(int) * 1);
	cudaMemset(d_counter_B, 0, sizeof(int));

	cudaMalloc(&d_adj_labels_B,      sizeof(int) * (m_B));
	cudaMalloc(&d_sec_adj_labels_B,      sizeof(int) * (m_B));
	cudaMalloc(&d_if_update_even_B,sizeof(bool) * n_A);
	cudaMalloc(&d_if_update_odd_B,sizeof(bool) * n_B);
	cudaMalloc(&d_r_offsets_B,        sizeof(E) * (n_B + 1));
	if(this->myconf.if_inverse==true){
		cudaMalloc(&d_neighborkey_out_B,sizeof(int) * (m_B));
		cudaMalloc(&d_neighborval_out_B,sizeof(int) * (m_B));
	}
	errorCheck("Out of Memory!");
	// copy to Global
	cudaSetDevice(0);
	cudaMemcpy(this->d_neighbors_A, &(neighborA[0]),sizeof(V)*m_A, cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_offsets_A, &(offsetA[0]), sizeof(E) * (n_A+1), cudaMemcpyHostToDevice);
	cudaSetDevice(1);
	cudaMemcpy(this->d_neighbors_B, &(neighborB[0]),sizeof(V)*(m_B), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_offsets_B, &(offsetB[0]), sizeof(E) * (n_B+1), cudaMemcpyHostToDevice);
}

template<typename V, typename E>
void MultiGprop<V, E>::postprocess()
{
	cudaSetDevice(0);
	cudaFree(d_neighbors_A);
	cudaFree(d_vertices_A);
	cudaFree(d_offsets_A);
	cudaFree(d_labels_A);
	cudaFree(d_sec_labels_even_A);
	cudaFree(d_sec_labels_odd_A);

	cudaFree(d_counter_A);
	cudaFree(d_if_update_even_A);
	cudaFree(d_if_update_odd_A);
	cudaFree(d_adj_labels_A);
	cudaFree(d_sec_adj_labels_A);
	if(this->myconf.if_inverse == true){
			cudaFree(d_neighborkey_out_A);
			cudaFree(d_neighborval_out_A);
	}

	cudaSetDevice(1);
	cudaFree(d_neighbors_B);
	cudaFree(d_vertices_B);
	cudaFree(d_offsets_B);
	cudaFree(d_labels_B);
	cudaFree(d_sec_labels_even_B);
	cudaFree(d_sec_labels_odd_B);

	cudaFree(d_counter_B);
	cudaFree(d_if_update_even_B);
	cudaFree(d_if_update_odd_B);
	cudaFree(d_adj_labels_B);
	cudaFree(d_sec_adj_labels_B);
	if(this->myconf.if_inverse == true){
		cudaFree(d_neighborkey_out_B);
		cudaFree(d_neighborval_out_B);
	}
}

template<typename V, typename E>
double MultiGprop<V,E>::run(int niter){
	double running_time;
	//query multigpu
	GPUquery();
	//assign data to GPU 0 and GPU 1
	preprocess();
	// lp
	running_time = perform_lp(G->n,G->m,niter);
	return running_time;
}

template<typename V, typename E>
double MultiGprop<V, E>::perform_lp(V n,E m,int niter)
{
	//init flag as false
	bool flag = false;
	// init odd_even_mask as true
	bool odd_even_mask = true;
	int changecount = 1000;// !=0
	const int VT = 7;
	const int VT_m = 1;
	const int nthreads = 128;
	const int n_blocks = divup(n,nthreads);
	const int m_blocks = divup(m,nthreads);

	const int n_blocks_A = divup(n_A,nthreads);
	const int n_blocks_B = divup(n_B,nthreads);
	const int m_blocks_A = divup(m_A,nthreads);
	const int m_blocks_B = divup(m_B,nthreads);
	const int m_blocks_div_A = divup(m_blocks_A,VT_m);
	const int m_blocks_div_B = divup(m_blocks_B,VT_m);

	const int nt_big = 32*4;
	int num_huge_A = 0;
	int num_big_A = 0;
	int num_medium_A = 0;
	int num_small_A = 0;
	int num_tiny_A = n_A;
	int num_huge_B = 0;
	int num_big_B = 0;
	int num_medium_B = 0;
	int num_small_B = 0;
	int num_tiny_B = n_B;
	// record the index in vectices array.
	int huge_index_A = 0;
	int big_index_A = 0;
	int medium_index_A = 0;
	int small_index_A = 0;
	int huge_index_B = 0;
	int big_index_B = 0;
	int medium_index_B = 0;
	int small_index_B = 0;


	const int layer = 4;
	const E medium = 32*layer;
	const E huge = 10000;
	cudaSetDevice(0);
	initialize_labels<<<n_blocks, nthreads>>>(d_labels_A, n);
	cudaSetDevice(1);
	initialize_labels<<<n_blocks, nthreads>>>(d_labels_B, n);
	Timer t0;
	t0.start();
	std::vector<V> vertices;
	vertices.resize(G->n);
	for (int i=0;i<G->n;++i) {
		vertices[i] = i;
	}
	// computer nodes size for both GPUs
	huge_index_A =  binary_search( 0, n_A-1, offsetA, vertices, huge);
	printf("A huge_index %d \n",huge_index_A);
	big_index_A =  binary_search( 0, n_A-1, offsetA, vertices, medium);
	printf("A big_index %d \n",big_index_A);
	medium_index_A = binary_search( 0, n_A-1, offsetA, vertices, 32);
	printf("A medium_index %d \n",medium_index_A);
	small_index_A = binary_search( 0, n_A-1, offsetA, vertices, 16);
	printf("A small_index %d \n",small_index_A);
	num_huge_A = huge_index_A+1;
	num_big_A = big_index_A-huge_index_A;
	num_medium_A = medium_index_A-big_index_A;
	num_small_A = small_index_A-medium_index_A;
	num_tiny_A = n_A-1-small_index_A;
	printf("GPU0, num_huge %d,num_big %d,num_medium %d,num_small %d,num_tiny %d \n",num_huge_A,num_big_A,num_medium_A,num_small_A,num_tiny_A);
	huge_index_B =  binary_search( 0, n_B-1, offsetB, vertices,  huge);
	printf("B huge_index %d \n",huge_index_B);
	big_index_B =  binary_search( 0,  n_B-1, offsetB, vertices, medium);
	printf("B big_index %d \n",big_index_B);
	medium_index_B = binary_search( 0, n_B-1, offsetB, vertices, 32);
	printf("B medium_index %d \n",medium_index_B);
	small_index_B = binary_search( 0,  n_B-1, offsetB, vertices, 16);
	printf("B small_index %d \n",small_index_B);
	num_huge_B = huge_index_B+1;
	num_big_B = big_index_B-huge_index_B;
	num_medium_B = medium_index_B-big_index_B;
	num_small_B = small_index_B-medium_index_B;
	num_tiny_B = n_B-1-small_index_B;
	printf("GPU1: num_huge %d,num_big %d,num_medium %d,num_small %d,num_tiny %d \n",num_huge_B,num_big_B,num_medium_B,num_small_B,num_tiny_B);
	//compute number of blocks for huge nodes
	const int nt_huge = 32*8;
	int nb_huge_A = 0;
	int nb_huge_B = 0;
	//huge nodes and some big nodes handled by GPU 0, others by GPU 1.
//	// if huge nide exists
	errorCheck("partition nodes");
	cudaSetDevice(0);
	if(num_huge_A>0){
		int* nb_huge_pn_A;
		nb_huge_pn_A = &nb_huge_A;
		cudaMalloc(&d_block_num_A,      sizeof(int) * (num_huge_A+1));
		//compute kernel
		compute_num_blocks<nt_huge><<<divup(num_huge_A,nthreads),nthreads>>>(d_offsets_A,0,huge_index_A, d_block_num_A);
		 //exclusive scan to compute nb_huge.
		thrust::exclusive_scan(thrust::device,d_block_num_A,d_block_num_A+num_huge_A+1,d_block_num_A);
		cudaMemcpy(nb_huge_pn_A,d_block_num_A+num_huge_A,sizeof(int), cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		// memory assignment after compute nb_huge.
		cudaMalloc(&d_block_v_A,      sizeof(int) * nb_huge_A);
		cudaMalloc(&d_block_id_A,      sizeof(int) * nb_huge_A);
		//assign each block its vertex
		assign_blocks<<<64, 128, 0, 0>>>(d_block_num_A, 0,huge_index_A, d_block_v_A,d_block_id_A);
		// global count for across block reducing
		cudaMalloc(&d_max_count_huge_A, sizeof(unsigned long long)*num_huge_A);
		cudaMemset(d_max_count_huge_A, 0, sizeof(unsigned long long)*num_huge_A);
		printf("huge node assignment success in GPU 0\n");
	}
	cudaSetDevice(1);
	if(num_huge_B>0){
		int* nb_huge_pn_B;
		nb_huge_pn_B = &nb_huge_B;
		cudaMalloc(&d_block_num_B,      sizeof(int) * (num_huge_B+1));
		//compute kernel
		compute_num_blocks<nt_huge><<<divup(num_huge_A,nthreads),nthreads>>>(d_offsets_B,0,huge_index_B, d_block_num_B);
		 //exclusive scan to compute nb_huge.
		thrust::exclusive_scan(thrust::device,d_block_num_B,d_block_num_B+num_huge_B+1,d_block_num_B);
		cudaMemcpy(nb_huge_pn_B,d_block_num_B+num_huge_B,sizeof(int), cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		// memory assignment after compute nb_huge.
		cudaMalloc(&d_block_v_B,      sizeof(int) * nb_huge_B);
		cudaMalloc(&d_block_id_B,      sizeof(int) * nb_huge_B);
		//assign each block its vertex
		assign_blocks<<<64, 128, 0, 0>>>(d_block_num_B, 0,huge_index_B, d_block_v_B,d_block_id_B);
		// global count for across block reducing
		cudaMalloc(&d_max_count_huge_B, sizeof(unsigned long long)*num_huge_B);
		cudaMemset(d_max_count_huge_B, 0, sizeof(unsigned long long)*num_huge_B);
		printf("huge node assignment success in GPU 1\n");
	}
	errorCheck("assign huge nodes");
	//allocate global memory for huge and big nodes
	if(num_big_A){
		cudaSetDevice(0);
		cudaMalloc(&(d_tables_A.keys),sizeof(int)*offsetA[big_index_A]);
		cudaMalloc(&(d_tables_A.vals),sizeof(int)*offsetA[big_index_A]);
	}
	if(num_big_B){
		cudaSetDevice(1);
		cudaMalloc(&(d_tables_B.keys),sizeof(int)*offsetB[big_index_B]);
		cudaMalloc(&(d_tables_B.vals),sizeof(int)*offsetB[big_index_B]);
	}
	errorCheck("hash table assignment");
//	// compute tiny nodes, how many edges for each warps
	int warpnumber_A = 1;
	int warpnumber_B = 1;
	std::vector<int> warp_v_A;
	std::vector<int> warp_begin_A;
	std::vector<int> warp_v_B;
	std::vector<int> warp_begin_B;
	int * d_warp_v_A;
	int * d_warp_begin_A;
	int * d_warp_v_B;
	int * d_warp_begin_B;
	int start_A = small_index_A;
	int start_B = small_index_B;
	if(small_index_A == -1){
		//all node are tiny
		start_A = 0;
	}
	warp_begin_A.push_back(start_A);
	for (int i = 0; i<= num_tiny_A;++i){
		if(offsetA[small_index_A+ i+1]-offsetA[start_A]>32){

			warpnumber_A++;
			for(int k =0;k<32-(offsetA[small_index_A+ i]-offsetA[start_A]);++k){
				warp_v_A.push_back(-1);
			}
			for(int k = 0; k<offsetA[small_index_A+ i+1] - offsetA[small_index_A+ i];++k){
				warp_v_A.push_back(small_index_A+ i);
			}
			start_A = small_index_A+ i;
			warp_begin_A.push_back(offsetA[start_A]);
		}else{

			for(int k = 0; k<offsetA[small_index_A+ i+1] - offsetA[small_index_A+ i];++k){
				warp_v_A.push_back(small_index_A+ i);
			}
		}
	}
	while(warpnumber_A*32-(int)warp_v_A.size()>0){
		warp_v_A.push_back(-1);
	}

	if(small_index_B == -1){
		//all node are tiny
		start_B = 0;
	}
	warp_begin_B.push_back(start_B);
	for (int i = 0; i<= num_tiny_B;++i){
		if(offsetB[small_index_B+ i+1]-offsetB[start_B]>32){

			warpnumber_B++;
			for(int k =0;k<32-(offsetA[small_index_B+ i]-offsetA[start_B]);++k){
				warp_v_B.push_back(-1);
			}
			for(int k = 0; k<offsetB[small_index_B+ i+1] - offsetA[small_index_B+ i];++k){
				warp_v_B.push_back(small_index_B+ i);
			}
			start_B = small_index_B+ i;
			warp_begin_B.push_back(offsetB[start_B]);
		}else{

			for(int k = 0; k<offsetB[small_index_B+ i+1] - offsetB[small_index_B+ i];++k){
				warp_v_B.push_back(small_index_B+ i);
			}
		}
	}
	while(warpnumber_B*32-(int)warp_v_B.size()>0){
		warp_v_B.push_back(-1);
	}

	int nt_tiny = 64;
	int nb_tiny2_A = divup(warpnumber_A, nt_tiny/32);
	int nb_tiny2_B = divup(warpnumber_B, nt_tiny/32);
	cudaSetDevice(0);
	cudaMalloc (&d_warp_v_A,sizeof(int)*(warpnumber_A*32));
	cudaMalloc (&d_warp_begin_A,sizeof(int)*warpnumber_A);
	cudaMemcpy(d_warp_v_A,&warp_v_A[0],sizeof(int)*(warpnumber_A*32),cudaMemcpyHostToDevice);
	cudaMemcpy(d_warp_begin_A,&warp_begin_A[0],sizeof(int)*warpnumber_A,cudaMemcpyHostToDevice);
	cudaSetDevice(1);
	cudaMalloc (&d_warp_v_B,sizeof(int)*(warpnumber_B*32));
	cudaMalloc (&d_warp_begin_B,sizeof(int)*warpnumber_B);
	cudaMemcpy(d_warp_v_B,&warp_v_B[0],sizeof(int)*(warpnumber_B*32),cudaMemcpyHostToDevice);
	cudaMemcpy(d_warp_begin_B,&warp_begin_B[0],sizeof(int)*warpnumber_B,cudaMemcpyHostToDevice);
	printf("tiny nodes assignment success \n");
	errorCheck("5st");
	const int nt_medium  = 32*4;
	// 1 node 1 warp
	const int nb_medium_A  = divup(num_medium_A, nt_medium/32*VT);
	const int nt_small = 32*2;
	const int nb_small_A = divup(num_small_A, nt_small/32*VT);
	const int nb_tiny_A  =divup(num_tiny_A, nt_tiny*VT);

	const int nb_medium_B  = divup(num_medium_B, nt_medium/32*VT);
	const int nb_small_B = divup(num_small_B, nt_small/32*VT);
	const int nb_tiny_B  =divup(num_tiny_B, nt_tiny*VT);


	cudaSetDevice(0);
	cudaMemcpy(this->d_vertices_A, &(vertices[0]), sizeof(V) * (G->n), cudaMemcpyHostToDevice);
	cudaSetDevice(1);
	cudaMemcpy(this->d_vertices_B, &(vertices[0]), sizeof(V) * (G->n), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	float t1(0.0);
	Timer t2;
	Timer t3;
	Timer t5;
	float t4(0.0);
	float t6(0.0);
	float t10(0.0);
	cudaDeviceSynchronize();
	t0.stop();
	std::cout << " precompute: " << t0.elapsed_time()<<std::endl;
	Timer t7;
	t7.start();
	cudaSetDevice(0);
	first_iteration<<<n_blocks_A, nthreads>>>(d_neighbors_A,d_offsets_A,d_labels_A, n_A);
	cudaSetDevice(1);
	first_iteration<<<n_blocks_B, nthreads>>>(d_neighbors_B,d_offsets_B,d_labels_B, n_B);
	cudaDeviceSynchronize();
	errorCheck("first iter");
	t7.stop();
	std::cout << "in interation: " << 0<< " running time: " << t7.elapsed_time() <<std::endl;
	bool onestep = myconf.onestep;
	bool if_inverse = myconf.if_inverse;
	int	changeiter = myconf.changeiter;
	bool opt_load = myconf.opt_load;
	Timer t9;
	    //the other iteration
	    for(int j = 1;j<niter;++j){
	    	// communication
	    	t9.start();
	    	if(flag==false){
	    		cudaMemcpy(d_labels_even_B,d_labels_even_A,sizeof(int)*n_A,cudaMemcpyDeviceToDevice);
				cudaMemcpy(d_labels_odd_A,d_labels_odd_B,sizeof(int)*n_B,cudaMemcpyDeviceToDevice);
	    	}else{
				if(onestep == true){
					// copy d_labels_even_A from GPU0 to GPU1
					// copy d_label_odd_B from GPU1 to GPU0
					cudaMemcpy(d_labels_even_B, d_labels_even_A, sizeof(int)*n_A, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_if_update_even_B, d_if_update_even_A, sizeof(bool)*n_A, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_labels_odd_A, d_labels_odd_B, sizeof(int)*n_B, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_if_update_odd_A, d_if_update_odd_B, sizeof(bool)*n_B, cudaMemcpyDeviceToDevice);
				}else{
					if(odd_even_mask){
						cudaMemcpy(d_labels_even_B,d_labels_even_A,sizeof(int)*n_A,cudaMemcpyDeviceToDevice);
						cudaMemcpy(d_labels_odd_A,d_labels_odd_B,sizeof(int)*n_B,cudaMemcpyDeviceToDevice);
						cudaMemcpy(d_if_update_even_B,d_if_update_even_A,sizeof(bool)*n_A,cudaMemcpyDeviceToDevice);
						cudaMemcpy(d_if_update_odd_A,d_if_update_odd_B,sizeof(bool)*n_B,cudaMemcpyDeviceToDevice);
					}else{
						cudaMemcpy(d_sec_labels_even_B,d_sec_labels_even_A,sizeof(int)*n_A,cudaMemcpyDeviceToDevice);
						cudaMemcpy(d_sec_labels_odd_A,d_sec_labels_odd_B,sizeof(int)*n_B,cudaMemcpyDeviceToDevice);
						cudaMemcpy(d_if_update_even_B,d_if_update_even_A,sizeof(bool)*n_A,cudaMemcpyDeviceToDevice);
						cudaMemcpy(d_if_update_odd_A,d_if_update_odd_B,sizeof(bool)*n_B,cudaMemcpyDeviceToDevice);
					}
				}
	    	}
	    	cudaDeviceSynchronize();
	    	errorCheck("copy memory");
	    	t9.stop();
	    	t2.start();
	    	t5.start();
	    	// label gathering
	    	if(flag==false){
	//    		gather_labels<int,int><<<m_blocks, nthreads, 0, 0>>>(d_neighbors, d_labels, d_adj_labels, m);
	    		cudaSetDevice(0);
	    		gather_labels_double<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighbors_A, d_labels_even_A, d_labels_odd_A, d_adj_labels_A, m_A);
	    		cudaSetDevice(1);
	    		gather_labels_double<int,int,VT_m><<<m_blocks_div_B, nthreads>>>(d_neighbors_B, d_labels_even_B, d_labels_odd_B, d_adj_labels_B, m_B);
	    	}else{
	    		if(onestep == true){
	    			if(!if_inverse){
	    				cudaSetDevice(0);
	    				gather_labels_double<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighbors_A, d_labels_even_A, d_labels_odd_A, d_adj_labels_A, m_A, d_if_update_even_A, d_if_update_odd_A);
	    				cudaSetDevice(1);
	    				gather_labels_double<int,int,VT_m><<<m_blocks_div_B, nthreads>>>(d_neighbors_B, d_labels_even_B, d_labels_odd_B, d_adj_labels_B, m_B, d_if_update_even_B, d_if_update_odd_B);
	    			}else{
	    				cudaSetDevice(0);
	    				gather_labels_double_inverse<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighborkey_out_A, d_neighborval_out_A, d_labels_even_A, d_labels_odd_A, d_adj_labels_A, m_A,d_if_update_even_A, d_if_update_odd_A);
	    				cudaSetDevice(1);
	    				gather_labels_double_inverse<int,int,VT_m><<<m_blocks_div_B, nthreads>>>(d_neighborkey_out_B, d_neighborval_out_B, d_labels_even_B, d_labels_odd_B, d_adj_labels_B, m_B,d_if_update_even_B, d_if_update_odd_B);
	    			}
	    		}else{
	    			if(changecount >0){
						if(odd_even_mask == true){
							// if true, generate d_adj_labels.
							if(!if_inverse){
								cudaSetDevice(0);
								gather_labels_double<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighbors_A, d_labels_even_A,d_labels_odd_A, d_adj_labels_A, m_A, d_if_update_even_A, d_if_update_odd_A);
								cudaSetDevice(1);
								gather_labels_double<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighbors_B, d_labels_even_B,d_labels_odd_B, d_adj_labels_B, m_B, d_if_update_even_B, d_if_update_odd_B);
							}else{
								cudaSetDevice(0);
								gather_labels_double_inverse<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighborkey_out_A,d_neighborval_out_A, d_labels_even_A, d_labels_odd_A, d_adj_labels_A, m_A,d_if_update_even_A,d_if_update_odd_A);
								cudaSetDevice(1);
								gather_labels_double_inverse<int,int,VT_m><<<m_blocks_div_B, nthreads>>>(d_neighborkey_out_B,d_neighborval_out_B, d_labels_even_B, d_labels_odd_B, d_adj_labels_B, m_B,d_if_update_even_B,d_if_update_odd_B);
							}
						}else{
							//if false, generate d_sec_adj_labels.
							if(!if_inverse){
								cudaSetDevice(0);
								gather_labels_double<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighbors_A, d_sec_labels_even_A,d_sec_labels_odd_A, d_sec_adj_labels_A, m_A, d_if_update_even_A, d_if_update_odd_A);
								cudaSetDevice(1);
								gather_labels_double<int,int,VT_m><<<m_blocks_div_B, nthreads>>>(d_neighbors_B, d_sec_labels_even_B,d_sec_labels_odd_B, d_sec_adj_labels_B, m_B, d_if_update_even_B, d_if_update_odd_B);
							}else{
								cudaSetDevice(0);
								gather_labels_double_inverse<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighborkey_out_A,d_neighborval_out_A, d_sec_labels_even_A,d_sec_labels_odd_A, d_sec_adj_labels_A, m_A, d_if_update_even_A, d_if_update_odd_A );
								cudaSetDevice(1);
								gather_labels_double_inverse<int,int,VT_m><<<m_blocks_div_A, nthreads>>>(d_neighborkey_out_B,d_neighborval_out_B, d_sec_labels_even_B,d_sec_labels_odd_B, d_sec_adj_labels_B, m_B, d_if_update_even_B, d_if_update_odd_B );
							}
						}
	    			}
	    		}
	    	}

	    	cudaDeviceSynchronize();
	    	t5.stop();
	    	errorCheck("loading label");

	    	if(flag == true){
	    		cudaSetDevice(0);
	    		cudaMemset(d_if_update_even_A, 0, n_A*sizeof(bool));
	    		cudaSetDevice(1);
	    		cudaMemset(d_if_update_odd_B, 0, n_B*sizeof(bool));
	    	}

	    	const int buffer = 1333;
			const int d = 3;
			const int w = 1333/d;
	     	cudaDeviceSynchronize();
	    	t3.start();
	    	// reset global hash table.
			if(num_big_A>0){
				cudaSetDevice(0);
				cudaMemset((d_tables_A.keys),0,sizeof(int)*offsetA[big_index_A]);
				cudaMemset((d_tables_A.vals),0,sizeof(int)*offsetA[big_index_A]);
				cudaSetDevice(1);
				cudaMemset((d_tables_B.keys),0,sizeof(int)*offsetB[big_index_B]);
				cudaMemset((d_tables_B.vals),0,sizeof(int)*offsetB[big_index_B]);
			}
	    	if(!opt_load||j<changeiter){
	    		//prepare for two step look up.
				if( opt_load&&onestep== false &&j== changeiter-1){
					cudaMemcpy(d_sec_labels_even_A, d_labels_even_A, sizeof(int)*n_A, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_sec_labels_odd_A, d_labels_odd_A, sizeof(int)*n_A, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_sec_labels_even_B, d_labels_even_B, sizeof(int)*n_B, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_sec_labels_odd_B, d_labels_odd_B, sizeof(int)*n_B, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_sec_adj_labels_A,d_adj_labels_A,sizeof(int)*m_A, cudaMemcpyDeviceToDevice);
					cudaMemcpy(d_sec_adj_labels_B,d_adj_labels_B,sizeof(int)*m_B, cudaMemcpyDeviceToDevice);
				}

				cudaDeviceSynchronize();
				cudaSetDevice(0);
				if(num_huge_A>0){
					 l_huge_update_syn2<32,nt_huge><<<nb_huge_A, nt_huge>>>( d_neighbors_A,d_offsets_A, d_adj_labels_A,d_vertices_A,d_block_v_A,d_block_id_A,d_max_count_huge_A,d_tables_A);
					 l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_A,nthreads),nthreads>>>(d_max_count_huge_A,  d_labels_even_A,d_counter_A, huge_index_A);
				}
				if(num_big_A>0){
					l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_A+VT-1)/VT, nt_big>>>(huge_index_A,big_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A, d_tables_A);
				}
				l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_A, nt_medium>>>(big_index_A,medium_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A,medium);
				l_small_update_syn<32,VT,nt_small><<<nb_small_A, nt_small>>>(medium_index_A,small_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A);
				l_tiny_update_syn2<<<nb_tiny2_A, nt_tiny,0, 0>>>(d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A, d_counter_A,d_warp_v_A,d_warp_begin_A,warpnumber_A);

				cudaSetDevice(1);
				if(num_huge_B>0){
					 l_huge_update_syn2<32,nt_huge><<<nb_huge_B, nt_huge>>>( d_neighbors_B,d_offsets_B, d_adj_labels_B,d_vertices_B,d_block_v_B,d_block_id_B,d_max_count_huge_B,d_tables_B);
					 l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_B,nthreads),nthreads>>>(d_max_count_huge_B,  d_labels_even_B,d_counter_B, huge_index_B);
				}
				if(num_big_B>0){
					l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_B+VT-1)/VT, nt_big>>>(huge_index_B,big_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B,d_vertices_B, d_counter_B, d_tables_B);
				}
				l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_B, nt_medium>>>(big_index_B,medium_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B,d_vertices_B, d_counter_B,medium);
				l_small_update_syn<32,VT,nt_small><<<nb_small_B, nt_small>>>(medium_index_B,small_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B,d_vertices_B, d_counter_B);
				l_tiny_update_syn2<<<nb_tiny2_B, nt_tiny,0, 0>>>(d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B, d_counter_B,d_warp_v_B,d_warp_begin_B,warpnumber_B);
	    	}else{
	    		flag = true;
	    		// look one step further
				if(onestep == true){
					cudaSetDevice(0);
					if(num_huge_A>0){
						l_huge_update_syn2<32,nt_huge><<<nb_huge_A, nt_huge>>>( d_neighbors_A,d_offsets_A, d_adj_labels_A,d_vertices_A,d_block_v_A,d_block_id_A,d_max_count_huge_A,d_tables_A);
						l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_A,nthreads),nthreads,0,0>>>(d_max_count_huge_A ,  d_labels_even_A ,d_counter_A , huge_index_A ,d_if_update_even_A);
					}
					if(num_big_A>0){
						l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_A+VT-1)/VT, nt_big>>>(huge_index_A, big_index_A ,d_neighbors_A, d_offsets_A, d_adj_labels_A, d_labels_even_A,d_vertices_A, d_counter_A, d_tables_A,d_if_update_even_A);
					}
					l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_A, nt_medium>>>(big_index_A,medium_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A,medium,d_if_update_even_A);
					l_small_update_syn<32,VT,nt_small><<<nb_small_A, nt_small>>>(medium_index_A,small_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A,d_if_update_even_A);
					l_tiny_update_syn2<<<nb_tiny2_A, nt_tiny,0, 0>>>(d_neighbors_A, d_offsets_A, d_adj_labels_A,d_labels_even_A, d_counter_A,d_warp_v_A,d_warp_begin_A,warpnumber_A,d_if_update_even_A);
					cudaSetDevice(1);
					if(num_huge_B>0){
						l_huge_update_syn2<32,nt_huge><<<nb_huge_B, nt_huge>>>( d_neighbors_B,d_offsets_B, d_adj_labels_B,d_vertices_B,d_block_v_B,d_block_id_B,d_max_count_huge_B,d_tables_B);
						l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_B,nthreads),nthreads>>>(d_max_count_huge_B ,  d_labels_B ,d_counter_B , huge_index_B ,d_if_update_odd_B);
					}
					if(num_big_B>0){
						l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_B+VT-1)/VT, nt_big>>>(huge_index_B, big_index_B ,d_neighbors_B, d_offsets_B, d_adj_labels_B, d_labels_even_B,d_vertices_B, d_counter_B, d_tables_B,d_if_update_odd_B);
					}
					l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_B, nt_medium>>>(big_index_B,medium_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B,d_vertices_B, d_counter_B,medium,d_if_update_odd_B);
					l_small_update_syn<32,VT,nt_small><<<nb_small_B, nt_small>>>(medium_index_B,small_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B,d_vertices_B, d_counter_B,d_if_update_odd_B);
					l_tiny_update_syn2<<<nb_tiny2_B, nt_tiny,0, 0>>>(d_neighbors_B, d_offsets_B, d_adj_labels_B,d_labels_even_B, d_counter_B,d_warp_v_B,d_warp_begin_B,warpnumber_B,d_if_update_odd_B);
	    		}else{
	    			//look two steps further
	    			if(odd_even_mask){
	    				// first iteration using d_adj_labels, output d_sec_labels
	    				cudaSetDevice(0);
	    				if(num_huge_A>0){
							l_huge_update_syn2<32,nt_huge><<<nb_huge_A, nt_huge>>>( d_neighbors_A,d_offsets_A, d_adj_labels_A, d_vertices_A, d_block_v_A, d_block_id_A, d_max_count_huge_A, d_tables_A);
							l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_A,nthreads),nthreads>>>(d_max_count_huge_A , d_sec_labels_even_A, d_counter_A, huge_index_A, d_if_update_even_A);
	    				}
	    				if(num_big_A>0){
	    					l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_A+VT-1)/VT, nt_big>>>(huge_index_A,big_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A, d_sec_labels_even_A, d_vertices_A, d_counter_A, d_tables_A, d_if_update_even_A);
	    				}
						l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_A, nt_medium>>>(big_index_A,medium_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A, d_sec_labels_even_A, d_vertices_A, d_counter_A, medium, d_if_update_even_A);
						l_small_update_syn<32,VT,nt_small><<<nb_small_A, nt_small>>>(medium_index_A,small_index_A,d_neighbors_A, d_offsets_A, d_adj_labels_A,d_sec_labels_even_A,d_vertices_A, d_counter_A, d_if_update_even_A);
						l_tiny_update_syn2<<<nb_tiny2_A, nt_tiny,0, 0>>>(d_neighbors_A, d_offsets_A, d_adj_labels_A, d_sec_labels_even_A, d_counter_A, d_warp_v_A, d_warp_begin_A, warpnumber_A, d_if_update_even_A);
						cudaSetDevice(1);
						if(num_huge_B>0){
							l_huge_update_syn2<32,nt_huge><<<nb_huge_B, nt_huge>>>( d_neighbors_B, d_offsets_B, d_adj_labels_B, d_vertices_B, d_block_v_B, d_block_id_B, d_max_count_huge_B, d_tables_B);
							l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_B,nthreads),nthreads>>>(d_max_count_huge_B ,  d_sec_labels_odd_B ,d_counter_B, huge_index_B,d_if_update_odd_B);
						}
						if(num_big_B>0){
							l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_B+VT-1)/VT, nt_big, 0,0>>>(huge_index_B,big_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_sec_labels_odd_B,d_vertices_B, d_counter_B, d_tables_B,d_if_update_odd_B);
						}
						l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_B, nt_medium,0, 0>>>(big_index_B,medium_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_sec_labels_odd_B,d_vertices_B, d_counter_B,medium,d_if_update_odd_B);
						l_small_update_syn<32,VT,nt_small><<<nb_small_B, nt_small,0, 0>>>(medium_index_B,small_index_B,d_neighbors_B, d_offsets_B, d_adj_labels_B,d_sec_labels_odd_B,d_vertices_B, d_counter_B,d_if_update_odd_B);
						l_tiny_update_syn2<<<nb_tiny2_B, nt_tiny,0, 0>>>(d_neighbors_B, d_offsets_B,  d_adj_labels_B,d_sec_labels_odd_B, d_counter_B,d_warp_v_B,d_warp_begin_B,warpnumber_B,d_if_update_odd_B);
	    			}else{
	    				// second iteration using d_sec_adj_labels, output d_labels
	    				cudaSetDevice(0);
	    				if(num_huge_A > 0){
	    					l_huge_update_syn2<32,nt_huge><<<nb_huge_A, nt_huge>>>( d_neighbors_A,d_offsets_A, d_sec_adj_labels_A,d_vertices_A,d_block_v_A,d_block_id_A,d_max_count_huge_A,d_tables_A);
	    					l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_A,nthreads),nthreads,0,0>>>(d_max_count_huge_A , d_labels_A, d_counter_A, huge_index_A, d_if_update_even_A);
	    				}
	    				if(num_big_A>0){
	    					l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_A+VT-1)/VT, nt_big>>>(huge_index_A,big_index_A,d_neighbors_A, d_offsets_A, d_sec_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A, d_tables_A,d_if_update_even_A);
	    				}
	    				l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_A, nt_medium>>>(big_index_A,medium_index_A,d_neighbors_A, d_offsets_A, d_sec_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A,medium,d_if_update_even_A);
						l_small_update_syn<32,VT,nt_small><<<nb_small_A, nt_small>>>(medium_index_A,small_index_A,d_neighbors_A, d_offsets_A, d_sec_adj_labels_A,d_labels_even_A,d_vertices_A, d_counter_A,d_if_update_even_A);
						l_tiny_update_syn2<<<nb_tiny2_A, nt_tiny,0, 0>>>(d_neighbors_A, d_offsets_A, d_sec_adj_labels_A,d_labels_A, d_counter_A,d_warp_v_A,d_warp_begin_A,warpnumber_A,d_if_update_even_A);
	    				cudaSetDevice(1);
						if(num_huge_B > 0){
							l_huge_update_syn2<32,nt_huge><<<nb_huge_B, nt_huge>>>( d_neighbors_B,d_offsets_B, d_sec_adj_labels_B,d_vertices_B,d_block_v_B,d_block_id_B,d_max_count_huge_B,d_tables_B);
							l_sub_huge_update_syn<32,nthreads><<<divup(num_huge_B,nthreads),nthreads,0,0>>>(d_max_count_huge_B ,  d_labels_odd_B ,d_counter_B, huge_index_B,d_if_update_odd_B);
						}
						l_big_update_syn2<32,VT, nt_big,buffer><<<(num_big_B+VT-1)/VT, nt_big>>>(huge_index_B,big_index_B,d_neighbors_B, d_offsets_B, d_sec_adj_labels_B, d_labels_odd_B, d_vertices_B, d_counter_B, d_tables_B,d_if_update_odd_B);
						l_medium_update_syn<32,VT,nt_medium,layer><<<nb_medium_B, nt_medium>>>(big_index_B,medium_index_B,d_neighbors_B, d_offsets_B, d_sec_adj_labels_B, d_labels_odd_B, d_vertices_B, d_counter_B,medium,d_if_update_odd_B);
						l_small_update_syn<32,VT,nt_small><<<nb_small_B, nt_small>>>(medium_index_B,small_index_B,d_neighbors_B, d_offsets_B, d_sec_adj_labels_B,d_labels_odd_B,d_vertices_B, d_counter_B,d_if_update_odd_B);
						l_tiny_update_syn2<<<nb_tiny2_B, nt_tiny>>>(d_neighbors_B, d_offsets_B, d_sec_adj_labels_B,d_labels_odd_B, d_counter_B,d_warp_v_B,d_warp_begin_B,warpnumber_B,d_if_update_odd_B);
	    			}
	    			odd_even_mask = !odd_even_mask;
	    		}
			}
	    	cudaDeviceSynchronize();
	    	t3.stop();
	    	t2.stop();
	    	//printf function
	    	errorCheck("computing");
	        t1+=t2.elapsed_time();
	        t4+=t3.elapsed_time();
	        t6+=t5.elapsed_time();
	        t10+=t9.elapsed_time();
	        cudaDeviceSynchronize();
	        changecount = get_count();
	        std::cout << "in interation: " <<j<< ", node update : "<<changecount<< " total time(t1): "<<t2.elapsed_time()<<", processing time(t2): "<<t3.elapsed_time()<<", load label time: "<<t5.elapsed_time()<<std::endl;

	        printf("\n");
	    }

	    std::cout << "Mix, all time: "<<t10 + t1 + t7.elapsed_time()<< " propogation time: "<< t4+t7.elapsed_time()<< " loading label: " << t6<<" memory transfer: " << t10<<std::endl;
	    return t10 + t1 + t7.elapsed_time();
}
