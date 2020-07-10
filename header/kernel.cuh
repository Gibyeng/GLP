// -*- coding: utf-8 -*-

#pragma once
#include <stdio.h>
#include "common/range.cuh"
#include "method/countmin.cuh"
#include "method/frequent.cuh"
#include "method/buffer.cuh"
#include "cub/cub/cub.cuh"
#include <cooperative_groups.h>
#include <curand_kernel.h>
#define FULL_MASK 0xffffffff
using namespace cooperative_groups;
namespace cg = cooperative_groups;
__inline__ __device__ int get_gid()
{
    return threadIdx.x + blockIdx.x * blockDim.x;
}

template<typename V, typename E>
__global__ void gather_labels(const V *neighbors, const V *labels, V *dst, E m, bool *if_update)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < m)
	{
		auto my_neighbor = neighbors[j];
		if(if_update[my_neighbor]){
			auto my_label = labels[my_neighbor];
			dst[j] =my_label;
		}

	}

}


template<typename V, typename E>
__global__ void push_gather_labels(const V *neighborkey_out ,const V *r_offsets, const V *labels, V *dst, E n, bool *if_update)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int v = blockIdx.x;
	if(v < n )
	{
		const int begin =r_offsets[v];
		const int end = r_offsets[v+1];

		for(int j = begin;j<end;j+= blockDim.x){
			int k = j+ threadIdx.x;
			if(k<end){
				auto my_label = labels[v];
				if(if_update[v]){
					dst[neighborkey_out[k]] =my_label;
				}
			}
		}
	}

}


//neighborkey_out: write to , neighborval_out: read from

template<typename V, typename E>
__global__ void gather_labels_inverse( V *neighborkey_out,V *neighborval_out, const V *labels, V *dst, E m, bool *if_update)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < m)
	{
		auto val =neighborval_out[j];
		if(if_update[val]){
			auto key = neighborkey_out[j];
			auto my_label = labels[val];
			dst[key] =my_label;
		}

	}
}

template<typename V, typename E,int VT>
__global__ void gather_labels_double_inverse(V* neighborkey_out, V* neighborval_out, V* labels_even,V* labels_odd, V* adj_labels,E m, bool* if_update_even, bool* if_update_odd)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < m)
	{
		auto val =neighborval_out[j];
		if(val %2==0){
			int i = val/2;
			if(if_update_even[i]){
				auto key = neighborkey_out[j];
				auto my_label = labels_even[i];
				adj_labels[key] =my_label;
			}
		}else{
			int i = (val-1)/2;
			if(if_update_odd[i]){
				auto key = neighborkey_out[j];
				auto my_label = labels_odd[i];
				adj_labels[key] =my_label;
			}
		}
	}
}


template<typename V, typename E,int VT>
__global__ void gather_labels(const V *neighbors, const V *labels, V *dst, E m, bool *if_update)
{
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];
			if(if_update[my_neighbor]){
				auto my_label = labels[my_neighbor];
				dst[j] =my_label;
			}
		}
	}
}

template<typename V, typename E,int VT>
__global__ void gather_labels_double(const V * neighbors, const V * labels_even, const V * labels_odd, V * adj_labels, E m, bool *if_update_even, bool *if_update_odd)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < m)
	{
		auto my_neighbor = neighbors[j];
		if(my_neighbor%2==0){
			int i = my_neighbor/2;
			auto my_label = labels_even[i];
			adj_labels[j] =my_label;

		}else{
			int i = (my_neighbor-1)/2;
			auto my_label = labels_odd[i];
			adj_labels[j] =my_label;
		}
	}

}

template<typename V, typename E,int VT>
__global__ void gather_labels_double(const V * neighbors, const V * labels_even, const V * labels_odd, V * adj_labels, E m)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < m)
	{
		auto my_neighbor = neighbors[j];
		if(my_neighbor%2==0){
			int i = my_neighbor/2;
			auto my_label = labels_even[i];
			adj_labels[j] =my_label;

		}else{
			int i = (my_neighbor-1)/2;
			auto my_label = labels_odd[i];
			adj_labels[j] =my_label;
		}
	}

}

template<typename V, typename E,int VT>
__global__ void gather_labels_cmp(const V *neighbors, const V *labels, V *dst, E m, bool *if_update,int* counter)
{
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];
			if(if_update[my_neighbor]){
				auto my_label = labels[my_neighbor];
				dst[j] =my_label;
				atomicAdd(counter,1);
			}

		}

	}

}
template<typename V, typename E,int VT>
__global__ void gather_labels_inverse( V *neighborkey_out,V *neighborval_out, const V *labels, V *dst, E m, bool *if_update)
{

	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto val =neighborval_out[j];
			if(if_update[val]){
				auto key = neighborkey_out[j];
				auto my_label = labels[val];
				dst[key] =my_label;
			}

		}

	}
}

template<typename V, typename E,int VT>
__global__ void gather_labels_inverse_cmp( V *neighborkey_out,V *neighborval_out, const V *labels, V *dst, E m, bool *if_update,int* counter )
{
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{

			auto val =neighborval_out[j];
			if(if_update[val]){
//				atomicAdd(counter,1);
//				auto key = neighborkey_out[j];
//				auto my_label = labels[val];
				dst[j%256] = val;
			}

		}

	}
}

// used
template<typename V, typename E>
__global__ void gather_labels(const V *neighbors, const V *labels, V *dst, E m)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x ;
	if(j < m)
	{
		auto my_neighbor = neighbors[j];

		auto my_label = labels[my_neighbor];

		dst[j] =my_label;

	}
}

//used
template<typename V, typename E, int VT>
__global__ void gather_labels(const V *neighbors, const V *labels, V *dst, E m)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];

			int my_label = labels[my_neighbor];
			dst[j] =my_label;
		}
	}

}
//one block one node, exact map
template<typename V, typename E,int NT>
__global__ void gather_labels(const V *neighbors,E* offsets, const V *labels, V *dst, V n, bool *if_update,int *map )
{

		int v= blockIdx.x  ;
		if(v < n && if_update[v]==true)
		{
			const int begin =  offsets[v];
			const int end =  offsets[v+1];
			for(int j = begin; j<end;j+= NT){
				if(j+threadIdx.x< end){
					dst[map[j+threadIdx.x]] = labels[v];
				}
			}
		}
}




template<typename V, typename E,int VT, int NT, int mapwidth>
__global__ void gather_labels(const V *neighbors, const V *labels, V *dst, E m, bool *if_update_bitmap)
{
	__shared__ bool s_map[mapwidth];
	// load to shared memory.
	for (int k = 0; k< mapwidth; k+=NT){
		int p = threadIdx.x+k;
		if(p<mapwidth){
			s_map[p] = if_update_bitmap[p];
		}

	}
	__syncthreads();
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		// rewrite back to global memory
		if(j < m)
		{
			auto my_neighbor = neighbors[j];
			if(s_map[hash(my_neighbor)%mapwidth]){
				auto my_label = labels[my_neighbor];
				dst[j] =my_label;
			}

		}
	}
}

// update some value in gcounterlist based on labellist and alpha.
template<typename V, typename E, int VT>
__global__ void ml_gather_counter(const V *labellist, V* gcounter,V* gcounterlist,E m,bool* gcounter_if_update)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto label = labellist[j];
			if(gcounter_if_update[label]){
				auto count = gcounter[label];
				gcounterlist[j] = count;
			}
		}
	}
}

// computer the difference of ..
template<typename V, typename E>
__global__ void compute_BL(V* d_gcounter_r,V* d_gcounter_w, int alpha ,E n, bool* gcounter_if_update)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < n)
	{
		int difference = d_gcounter_r[j] - d_gcounter_w[j];
		if(difference < -alpha || difference > alpha){
			gcounter_if_update[j]  = 1;
		}
	}
}

template<typename V, typename E,int NT>
__global__ void compute_candidate(const V *neighbors, E *offsets,bool *if_update, bool *candidate)
{
	const int v = blockIdx.x;
	const int begin = offsets[v];
	const int end = offsets[v + 1];

	if( if_update[v]==true){
		for (int i = begin; i < end; i += (int)NT) {
			if(i+(int)threadIdx.x< end){
				candidate[neighbors[i+(int)threadIdx.x]] = true;
			}
		}
	}

}
// one thread one node.
template<typename V, typename E>
__global__ void compute_map(const V *neighbors, E *offsets, int m,int *map,int *mapcount)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i<m){
		int key = neighbors[i];
		int count = mapcount[key];
		int begin = offsets[key];
		int end   = offsets[key+1];
		auto old = atomicCAS(&map[begin + count],0,i);
		while (old != 0 && begin + count < end ){
			count++;
			old = atomicCAS(&map[begin + count],0,i);

		}
		atomicAdd(&mapcount[key],1);
	}
}
//compute if a vertice 's label changed, and compute counter
__global__ void compute_if_update (int* label,int* label_temp,bool* if_update,int *counter,int n)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < n)
	{
		if(label[j]!=label_temp[j]){
			if_update[j] = true;
			atomicAdd(counter,1);
		}
	}
}
//compute how mang label changed
__global__ void compute_counter (int* label,int* label_temp,int *counter,int n)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	if(j < n)
	{
		if(label[j]!=label_temp[j]){
			atomicAdd(counter,1);
		}
	}
}
template<typename V, typename E,int NT>
__global__ void compare(bool *if_update, bool *if_update2)
{
	const int v = blockIdx.x;
	if(threadIdx.x ==0){
		if( if_update[v]!= if_update2[v]){
			printf(" v %d ",v);
		}
	}

}

__device__ int hash(int k) {
		   k ^= k >> 16;
		   k *= 0x85ebca6b;
		   k ^= k >> 13;
		   k *= 0xc2b2ae35;
		   k ^= k >> 16;
		   return abs(k);
}

__device__ int hash1(int k) {
        k ^= k >> 10;
        k *= 0x82ecba63;
        k ^= k >> 12;
        k *= 0xc3b15e32;
        k ^= k >> 9;
        return abs(k);
    }

__device__ int hash2(int k) {
        k ^= k >> 12;
        k *= 0xe85cb46c;
        k ^= k >> 17;
        k *= 0xc3b2cf55;
        k ^= k >> 12;
        return abs(k);
    }

__device__ int hash3(int k) {
        k ^= k >> 16;
        k *= 0x85ebca6b;
        k ^= k >> 13;
        k *= 0xc2b2ae35;
        k ^= k >> 16;
        return abs(k);
    }

template<typename V>
__global__ void initialize_labels(V *labels, V n)
{
    for (auto i: grid_stride_range(n)) {
        // Vertex i is labeled i initially
        labels[i] = i;
    }
}

template<typename V>
__global__ void initialize_labels_by_seeds(V *labels, V n, bool* seeds)
{
    for (auto i: grid_stride_range(n)) {
        // seeds vertice i is labeled with i initially, non-rooted labeled with -1;
    	if(seeds[i] == false){
    		 labels[i] = -1;

    	}else{
    		 labels[i] = i;
    	}

    }
}

template<typename V>
__global__ void initialize_neighbor(V *neighbor_in, V m)
{
    for (auto i: grid_stride_range(m)) {
        // Vertex i is labeled i initially
        neighbor_in[i] = i;
    }
}



__global__ void KernelVersionShim(){}

// Lock-free hash tables on global memory
// - Just perform counting
// - Block per vertex
__global__ void count_lockfree
(int *adj_labels, int *offsets, GlobalHT g_tables)
{
    const int v = blockIdx.x;
    const int begin = offsets[v];
    const int end = offsets[v + 1];
    for (auto i: block_stride_range(begin, end)) {
        auto key = adj_labels[i] + 1;
        g_tables.increment(begin, end, key);
    }
}

template<typename V>
__global__ void gen_prob
(double* randprob, V n){
	for (auto i: grid_stride_range(n)) {
			curandState state;
			curand_init(clock64(), i, 0, &state);
			randprob[i] = curand_uniform(&state);
	    }
}

template<typename V>
__global__ void gen_prob2
(double* randprob, V n){
	for (auto i: grid_stride_range(n)) {
			randprob[i] = hash(i)/INT_MAX;
	    }
}

__global__ void update_labels
(uint32_t *keys, int *label_index, int n, int *labels, int *counter)
{
    int gid = get_gid();
    __shared__ int s_count;
    if (threadIdx.x == 0) {
        s_count = 0;
    }
    __syncthreads();

    if (gid < n) {
        int label = keys[label_index[gid]] - 1;
        if (label != labels[gid]) {
            atomicAdd(&s_count, 1);
        }
        labels[gid] = label;
    }

    if (threadIdx.x == 0) {
        atomicAdd(counter, s_count);
    }
}


template<int SC>
__device__ void flush_s2g
(GlobalHT &g_tables, int begin, int end,
 SharedHT<SC> &s_table, int &my_max_key, int &my_max_count)
{
    // Flush s_table to g_tables
    for (auto i: block_stride_range(SC)) {
        auto key = s_table.keys[i];
        auto count = s_table.vals[i];
        if (key > 0) {
            auto c = g_tables.increment(begin, end, key, count);
            if (c > my_max_count) {
                my_max_key = key;
                my_max_count = c;
            }
        }
    }
    __syncthreads();
    s_table.clear();
    __syncthreads();
}


// Kernel fusion and shared memory hash table
// - TS elements are first aggregated on the shared memory,
//   and then flush it to the global hash tables
// - Still block per vertex
template<int NT, int TS>
__global__ void ht_update
(int *neighbors,int *offsets, int *labels, GlobalHT g_tables, int *counter, int v_offset=0)
{
    constexpr int SC = TS + NT;  // Capacity of the smem hash table
    __shared__ SharedHT<SC> s_table;
    s_table.clear();
    __syncthreads();

    __shared__ int s_max_keys[NT];
    __shared__ int s_max_counts[NT];

    const int v = blockIdx.x;
    const int begin = offsets[v];
    const int end = offsets[v + 1];
    //if there is single node, key = my_max_key.
    int my_max_key = v+1;
    int my_max_count = 0;
    // one block start with begin end with end. NT thread per block.
    for (int i = begin; i < end; i += NT) {
        int j = i + threadIdx.x;
        if (j < end) {
            auto key = labels[neighbors[j]] + 1;
            s_table.increment(key);
        }
        // s_table is full.
        if (s_table.nitems >= SC - NT) {
            __syncthreads();
            flush_s2g(g_tables, begin, end, s_table, my_max_key, my_max_count);
        }
    }
    __syncthreads();
    if (s_table.nitems > 0) {
        flush_s2g(g_tables, begin, end, s_table, my_max_key, my_max_count);
    }
    s_max_keys[threadIdx.x] = my_max_key;
    s_max_counts[threadIdx.x] = my_max_count;
    __syncthreads();

    if (threadIdx.x == 0) {
        for (int i = 1; i < NT; ++i) {
            if (s_max_counts[i] > my_max_count) {
                my_max_key = s_max_keys[i];
                my_max_count = s_max_counts[i];
            }
        }
        auto lbl = labels[v + v_offset];
        if (lbl != my_max_key - 1) {
            atomicAdd(counter, 1);
            labels[v + v_offset] = my_max_key - 1;
        }
    }
}


//base line
template<int NT>
__global__ void update_lockfree
(int *neighbors, int *offsets, int *labels,int *temp_labels, GlobalHT g_tables, int* vertices,int *counter, double scale_factor=1.0)
{
    __shared__ int s_max_keys[NT];
    __shared__ int s_max_counts[NT];
    //const int v = blockIdx.x;
    const int v = vertices[blockIdx.x];
    const int begin = offsets[v];
    const int t_beg = begin * scale_factor;
    const int end = offsets[v + 1];
    const int t_end = end * scale_factor;
    int my_max_count = 0;
    int my_max_key = 0;
    for (auto i: block_stride_range(begin, end)) {
        auto label= labels[i];

        // Keys are >= 1; labels are >= 0
        auto key = label + 1;
        auto c = g_tables.increment(t_beg, t_end, key);
        if (c > my_max_count) {
            my_max_key = key;
            my_max_count = c;
        }
    }
    s_max_keys[threadIdx.x] = my_max_key;
    s_max_counts[threadIdx.x] = my_max_count;
    __syncthreads();

    if (threadIdx.x == 0) {
        for (int i = 1; i < NT; ++i) {  // Naive reduction
            if (s_max_counts[i] > my_max_count) {
                my_max_key = s_max_keys[i];
                my_max_count = s_max_counts[i];
            }
        }

        if (temp_labels[v] != my_max_key - 1) {
        	atomicAdd(counter, 1);
        	temp_labels[v] = my_max_key - 1;
        }

    }
}
//base line
template<int NT>
__global__ void update_lockfree
(int offset,int n,int *neighbors, int *offsets, int *labels,int *temp_labels, GlobalHT g_tables, int* vertices,int *counter, double scale_factor=1.0)
{
    __shared__ int s_max_keys[NT];
    __shared__ int s_max_counts[NT];

    if(blockIdx.x+ offset < n){
    	//const int v = blockIdx.x + offset;
		const int v = vertices[blockIdx.x+ offset];
		const int begin = offsets[v];
		const int t_beg = begin * scale_factor;
		const int end = offsets[v + 1];
		const int t_end = end * scale_factor;
		int my_max_count = 0;
		int my_max_key = 0;
		for (auto i: block_stride_range(begin, end)) {
			auto label= labels[i];

			// Keys are >= 1; labels are >= 0
			auto key = label + 1;
			auto c = g_tables.increment(t_beg, t_end, key);
			if (c > my_max_count) {
				my_max_key = key;
				my_max_count = c;
			}
		}
		s_max_keys[threadIdx.x] = my_max_key;
		s_max_counts[threadIdx.x] = my_max_count;
		__syncthreads();

		if (threadIdx.x == 0) {
			for (int i = 1; i < NT; ++i) {  // Naive reduction
				if (s_max_counts[i] > my_max_count) {
					my_max_key = s_max_keys[i];
					my_max_count = s_max_counts[i];
				}
			}

			if (temp_labels[v] != my_max_key - 1) {
				atomicAdd(counter, 1);
				temp_labels[v] = my_max_key - 1;
			}
		}
    }
}

// random pick a neighbor in its neighborhoods.
__global__ void first_iteration (int *neighbors,int *offsets,int *temp_labels, int n){
	const int v= get_gid();
	if(v < n){
		const int begin = offsets[v];
		const int end = offsets[v+1];

		if(begin <end){
			// random pick one neighbor
			int pick = hash(n)%(begin-end)+ begin;
			//In the first iteration, label[i] == i;
			temp_labels[v] = neighbors[pick];
		}else{
			temp_labels[v] = v;
		}

	}
}

// assign label by seed nodes,  push label from seed node.
__global__ void first_iteration (int *neighbors,int *offsets,int *temp_labels, int n,bool* seednodes){
	const int v= get_gid();
	if(v < n ){
		if(seednodes[v] == true){
			//seed node label =vid;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//
			for (int u_index = begin; u_index<end; ++u_index){
				int u = neighbors[offsets[u_index]];
//				temp_labels[u] = v;
			}
		}

	}
}

template<int TW,int NT>
__global__ void l_huge_update_syn
(int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int* block_v,int* block_id,int* max_count ,int *counter,SharedCM<200,3>* cm)
{
	__shared__ int s_max_key;
	__shared__ int s_max_count;
	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = 0;
	int inserted = 0;

	if (threadIdx.x == 0) {
		s_max_key = 0;
		s_max_count = 0;
	}
	__syncthreads();
	int j = begin+block_index*NT+threadIdx.x;
	if (j<end) {
		int key;
		if (j < end) {
			key = adj_labels[j]+1;

		}else{
			key = 0;
		}
		auto lmask = __match_any_sync(FULL_MASK,key);
		auto lleader  = __ffs(lmask) - 1;
		auto count  = __popc(lmask);
		if (g.thread_rank()== lleader&&key!=0){
			int mincount = cm[v].increment(key,count);
			int c =mincount;
			if(c >my_max_count){
				my_max_count = c;
				my_max_key   = key;
			}
		}
		//block reduce
		auto m = atomicMax(&s_max_count, my_max_count);
		if (m < my_max_count && s_max_count == my_max_count) {
			s_max_key = my_max_key;
		}
		__syncthreads();

		// Try to update the label of the vertex a block handles
		if (threadIdx.x == 0) {
			auto ret = atomicMax(&max_count[v], s_max_count);
			if (ret < s_max_count && max_count[v] == s_max_count) {
				auto lbl = temp_labels[v];
				temp_labels[v] = s_max_key - 1;

				if (lbl != s_max_key - 1) {
					// May count redundantly
//					atomicAdd(counter, 1);
				}
			}
		}
	}
}

// find the most frequent label
template<int TW,int NT>
__global__ void l_huge_update_syn2
(int *neighbors,int *offsets, int *adj_labels,int* vertices,int* block_v,int* block_id,unsigned long long* max_count ,GlobalHT gtable)
{
	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = 0;

	int j = begin+block_index*NT+threadIdx.x;
	int key;
	if (j < end) {
		key = adj_labels[j]+1;

	}else{
		key = 0;
	}
	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	if (g.thread_rank()== lleader&&key!=0){
		count = gtable.increment(begin,end,key,count);
		if(count >my_max_count){
			my_max_count = count;
			my_max_key   = key;
		}
	}

	//block reduce
	typedef cub::BlockReduce<unsigned long long,NT> BlockReduce;
	__shared__ typename BlockReduce::TempStorage temp_storage;
	unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;

	unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());
	__syncthreads();
	// Try to update the label of the vertex a block handles
	if (j<end) {
		if (threadIdx.x == 0 ) {
			 atomicMax(&max_count[v], aggregate);
		}
	}
}

template<int TW,int NT>
__global__ void l_huge_update_syn2_labelload
(int *neighbors,int *offsets, int *labels_read,int* vertices,int* block_v,int* block_id,unsigned long long* max_count ,GlobalHT gtable)
{
	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = 0;

	int j = begin+block_index*NT+threadIdx.x;
	int key;
	if (j < end) {
		key = labels_read[neighbors[j]]+1;

	}else{
		key = 0;
	}
	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	if (g.thread_rank()== lleader&&key!=0){
		count = gtable.increment(begin,end,key,count);

		if(count >my_max_count){
			my_max_count = count;
			my_max_key   = key;
		}
	}

	//block reduce
	typedef cub::BlockReduce<unsigned long long,NT> BlockReduce;
	__shared__ typename BlockReduce::TempStorage temp_storage;
	unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;

	unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());
	__syncthreads();
	// Try to update the label of the vertex a block handles
	if (j<end) {
		if (threadIdx.x == 0  ) {
			 atomicMax(&max_count[v], aggregate);
		}
	}
}

// for seed label prop.
template<int TW,int NT>
__global__ void l_huge_update_syn2_seed
(int *neighbors,int *offsets, int *adj_labels,int* vertices,int* block_v,int* block_id,unsigned long long* max_count ,GlobalHT gtable, bool* seednodes)
{
	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = 0;

	int j = begin+block_index*NT+threadIdx.x;
	int key;
	if (j < end) {
		key = adj_labels[j]+1;

	}else{
		key = -1;
	}
	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	bool seed = seednodes[v];

	if (g.thread_rank()== lleader&&key!=-1 && seed == false){
		count = gtable.increment(begin,end,key,count);
		if(count >my_max_count){
			my_max_count = count;
			my_max_key   = key;
		}
	}
	if(seed == true){
		my_max_key = adj_labels[v];
		my_max_count = 1;
		unsigned long long valkey =my_max_count <<32 + my_max_key;
		max_count[v] = valkey;
	}else{
		//block reduce
		typedef cub::BlockReduce<unsigned long long,NT> BlockReduce;
		__shared__ typename BlockReduce::TempStorage temp_storage;
		unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;

		unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());
		__syncthreads();
		// Try to update the label of the vertex a block handles
		if (j<end) {
			if (threadIdx.x == 0 ) {
				 atomicMax(&max_count[v], aggregate);
			}
		}
	}
}

// the score function
__device__ int compute_score(int labelcount, int v, int globalcount){
 const int theta = 2;
 return (theta+1)*labelcount - globalcount;
}

template<int TW,int NT>
__global__ void l_huge_update_syn2_seed_test
(int *neighbors,int *offsets, int *adj_labels,int* vertices,int* block_v,int* block_id,unsigned long long* max_count ,GlobalHT gtable, bool* seednodes)
{
	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_score = INT_MIN;

	int j = begin+block_index*NT+threadIdx.x;
	int key;
	if (j < end) {
		key = adj_labels[j]+1;

	}else{
		key = -1;
	}
	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	bool seed = seednodes[v];
	int score = INT_MIN;
	if (g.thread_rank()== lleader&&key!=-1 && seed == false){
		count = gtable.increment(begin,end,key,count);
		score = compute_score(count,v,0);
		if(score >my_max_score){
			my_max_score = score;
			my_max_key   = key;
		}
	}
	if(seed == true){
		my_max_key = adj_labels[v];
		my_max_score = 1;
		unsigned long long valkey =my_max_score <<32 + my_max_key;
		max_count[v] = valkey;
	}else{
		//block reduce
		typedef cub::BlockReduce<unsigned long long,NT> BlockReduce;
		__shared__ typename BlockReduce::TempStorage temp_storage;
		unsigned long long valkey = ((unsigned long long)my_max_score <<32 )+ my_max_key;

		unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());
		__syncthreads();
		// Try to update the label of the vertex a block handles
		if (j<end) {
			if (threadIdx.x == 0 ) {
				 atomicMax(&max_count[v], aggregate);
			}
		}
	}
}
// update label and counter
template<int TW,int NT>
__global__ void l_sub_huge_update_syn
(unsigned long long * max_count , int* temp_labels ,int *counter, int n)
{
	const int v = threadIdx.x + blockIdx.x * blockDim.x;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = 1;
	int lbl = 0;
	if(v < n){
		// cut the second part of array to get the label of v
		my_max_key = max_count[v];
		lbl = temp_labels[v];
	}
	auto mask = __ballot_sync(FULL_MASK,lbl != my_max_key - 1);
	auto c  = __popc(mask);
	if(g.thread_rank()==0){
//		atomicAdd(counter, c);
	}
	// update
	if (lbl != my_max_key - 1) {
		temp_labels[v] = my_max_key - 1;
	}
}

// mask v as 1 if v should be updated
template<int TW,int NT>
__global__ void l_sub_huge_update_syn
(unsigned long long * max_count , int *temp_labels ,int *counter, int n, bool * if_update)
{
	const int v = threadIdx.x + blockIdx.x * blockDim.x;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = 1;
	int lbl = 0;
	if(v < n){
		// cut the second part of array to get the label of v
		my_max_key = max_count[v];
		lbl = temp_labels[v];
	}
	auto mask = __ballot_sync(FULL_MASK,lbl != my_max_key - 1);
	auto c  = __popc(mask);
	if(g.thread_rank()==0){
//		atomicAdd(counter, c);
	}
	// update
	if (lbl != my_max_key - 1) {
		temp_labels[v] = my_max_key - 1;
		if_update[v] = 1;
	}
}


template<int TW,int NT>
__global__ void l_huge_update_syn_p
(int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int* block_v,int* block_id,int* max_count ,int *counter,SharedCM<200,3>* cm, double *prob,double *rand_prob)
{
	 __shared__ int s_max_key;
	 __shared__ int s_max_count;
	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = 0;
	int inserted = 0;

	if (threadIdx.x == 0) {
		s_max_key = 0;
		s_max_count = 0;
	}
	double p = 1-prob[v];
	if(rand_prob[v]<p){
		return;
	}
	__syncthreads();
	int j = begin+block_index*NT+threadIdx.x;
	if (j<end) {
		int key;
		if (j < end) {
			key = adj_labels[j]+1;

		}else{
			key = 0;
		}
		auto lmask = __match_any_sync(FULL_MASK,key);
		auto lleader  = __ffs(lmask) - 1;
		auto count  = __popc(lmask);
		if (g.thread_rank()== lleader&&key!=0){
			int mincount = cm[v].increment(key,count);
			int c =mincount;
			if(c >my_max_count){
				my_max_count = c;
				my_max_key   = key;
			}
		}
		//block reduce
		auto m = atomicMax(&s_max_count, my_max_count);
		if (m < my_max_count && s_max_count == my_max_count) {
			s_max_key = my_max_key;
		}
		__syncthreads();
		#pragma unroll
		// Try to update the label of the vertex a block handles
		if (threadIdx.x == 0) {
			auto ret = atomicMax(&max_count[v], s_max_count);
			if (ret < s_max_count && max_count[v] == s_max_count) {
				auto lbl = temp_labels[v];
				temp_labels[v] = s_max_key - 1;

				if (lbl != s_max_key - 1) {
					// May count redundantly
//					atomicAdd(counter, 1);
					// update
					prob[v] = 0;
				}else{
					prob[v] = 0.5+prob[v]/2;
				}
			}
		}
	}
}

// score = delta * local - global
// basic version
template<int TW,int NT>
__global__ void ml_huge_update_syn
(int *neighbors,int *offsets, int *adj_labels,int* vertices,int* block_v,int* block_id, long long* max_count ,GlobalHT gtable,int* gcounter,int delta)
{

	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = INT_MIN;
	int inserted = 0;


	__syncthreads();
	int j = begin+block_index*NT+threadIdx.x;

	int key;
	if (j < end) {
		key = adj_labels[j]+1;

	}else{
		key = 0;
	}
	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	int score = INT_MIN;

	if (g.thread_rank()== lleader&&key!=0){
		count = gtable.increment(begin,end,key,count);
		score = delta *count - gcounter[key-1];
		if(score >my_max_count){
			my_max_count = score;
			my_max_key   = key;
		}
	}

	//block reduce
	typedef cub::BlockReduce<long long,NT> BlockReduce;
	__shared__ typename BlockReduce::TempStorage temp_storage;
	long long valkey = ((long long)my_max_count <<32 )+ my_max_key;
	long long aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

	__syncthreads();
	#pragma unroll
	// Try to update the label of the vertex a block handles
	if (threadIdx.x == 0) {
		auto ret = atomicMax(&max_count[v], aggregate);
	}

}

// update label and counter
template<int TW,int NT>
__global__ void ml_sub_huge_update_syn
( long long * max_count , int* temp_labels ,int *counter, int n, int * gcounter_temp)
{
	const int v = threadIdx.x + blockIdx.x * blockDim.x;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = 1;
	int lbl = 0;
	if(v < n){
		// cut the second part of array to get the label of v
		my_max_key = max_count[v];
		lbl = temp_labels[v];
	}
	auto mask = __ballot_sync(FULL_MASK,lbl != my_max_key - 1);
	auto c  = __popc(mask);
	if(v < n){
		if(g.thread_rank()==0){
			atomicAdd(counter, c);
		}
		// update
		if (lbl != my_max_key - 1) {
			temp_labels[v] = my_max_key - 1;
		}
		atomicAdd(&gcounter_temp[my_max_key-1],1);
	}
}

// mask v as 1 if v should be updated
template<int TW,int NT>
__global__ void ml_sub_huge_update_syn
( long long * max_count , int* temp_labels ,int *counter, int n, int * gcounter_temp, bool* if_update)
{
	const int v = threadIdx.x + blockIdx.x * blockDim.x;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = 1;
	int lbl = 0;
	if(v < n){
		// cut the second part of array to get the label of v
		my_max_key = max_count[v];
		lbl = temp_labels[v];
	}
	auto mask = __ballot_sync(FULL_MASK,lbl != my_max_key - 1);
	auto c  = __popc(mask);
	if(v < n){
		if(g.thread_rank()==0){
			atomicAdd(counter, c);
		}
		// update
		if (lbl != my_max_key - 1) {
			temp_labels[v] = my_max_key - 1;
			if_update[v] = 1;
		}
		atomicAdd(&gcounter_temp[my_max_key-1],1);
	}
}

template<int TW,int NT>
__global__ void ml_huge_update_syn_labelload
(int *neighbors,int *offsets, int *labels,int* vertices,int* block_v,int* block_id, long long* max_count ,GlobalHT gtable,int* gcounter,int delta)
{

	const int v = block_v[blockIdx.x];
	const int begin = offsets[v];
	const int end = offsets[v+1];
	const int block_index = block_id[blockIdx.x];

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	int my_max_key = v+1;
	int my_max_count = INT_MIN;
	int inserted = 0;


	__syncthreads();
	int j = begin+block_index*NT+threadIdx.x;

	int key;
	if (j < end) {
		key = labels[neighbors[j]]+1;

	}else{
		key = 0;
	}
	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	int score = INT_MIN;

	if (g.thread_rank()== lleader&&key!=0){
		count = gtable.increment(begin,end,key,count);
		score = delta *count - gcounter[key-1];
		if(score >my_max_count){
			my_max_count = score;
			my_max_key   = key;
		}
	}

	//block reduce
	typedef cub::BlockReduce<long long,NT> BlockReduce;
	__shared__ typename BlockReduce::TempStorage temp_storage;
	long long valkey = ((long long)my_max_count <<32 )+ my_max_key;
	long long aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

	__syncthreads();
	#pragma unroll
	// Try to update the label of the vertex a block handles
	if (threadIdx.x == 0) {
		auto ret = atomicMax(&max_count[v], aggregate);
	}

}



template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter)
{

	 __shared__ int h1vals[W];
	 __shared__ int h2vals[W];
	 __shared__ int h3vals[W];


    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);

	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){

    				h1vals[threadIdx.x+stride] = 0;
					h2vals[threadIdx.x+stride] = 0;
					h3vals[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
			//const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key;
				if (j < end) {
					key = adj_labels[j]+1;
				}else{
					key = 0;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);

				if (g.thread_rank()== lleader&&key!=0){
					int h1 = hash1(key)%W;
					int h2 = hash2(key)%W;
					int h3 = hash3(key)%W;
					int mincount;
					int count1 = atomicAdd(&h1vals[h1],count)+count;
					int count2 = atomicAdd(&h2vals[h2],count)+count;
					int count3 = atomicAdd(&h3vals[h3],count)+count;
					mincount = min(count1,count2);
					mincount = min(mincount,count3);

					if(mincount >my_max_count){
						my_max_count = mincount;
						my_max_key   = key;
					}
				}
			}

			__syncthreads();
			typedef cub::BlockReduce<unsigned long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];

				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}


// W > distinct label number.
template<int TW, int VT,int NT,int htW,int W,int d>
__global__ void l_big_update_syn3
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable)
{
	__shared__  int vals[htW];
    __shared__  int keys[htW];
 	__shared__ SharedCM<W,d> s_cm;
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
     	s_cm.clear();
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
			#pragma unroll
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();

    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int cm_max_key = v+1;
			int my_max_count = 0;
			int cm_max_count = 0;
			bool inserted_HT = false;
			//first insert to shared memory
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);

				// hash to shared memory.
				auto hashkey = hash(key);
				hashkey = (hashkey)%W;
				auto old = keys[hashkey];
				if(old == 0 || old == key){
					old = atomicCAS(&keys[hashkey],0,key);
					if(old ==0|| old == key){
						//successful insert into shared memory
						inserted_HT = true;
						int tempcount = atomicAdd(&vals[hashkey],count);
						count += tempcount;
						if(my_max_count < count){
							my_max_count = count;
							my_max_key = key;
						}
					}
				}
				//hash to CM sketch
				if(inserted_HT == false){
					auto cm_count = s_cm.increment(key, count);
					if(cm_max_count < cm_count){
						cm_max_count = cm_count;
						cm_max_key = key;
					}
				}
			}
			__syncthreads();
			// find the most frequent label in CM sketch
			typedef cub::BlockReduce<unsigned long long,NT> BlockReduce4;
			__shared__ typename BlockReduce4::TempStorage temp_storage4;
			unsigned long long valkey = ((unsigned long long)cm_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce4(temp_storage4).Reduce(valkey,cub::Max());

			//reduce on maximum count of shared hash table
			typedef cub::BlockReduce<unsigned long long,NT> BlockReduce2;
			__shared__ typename BlockReduce2::TempStorage temp_storage2;
			unsigned long long valkey2 = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate2 = BlockReduce2(temp_storage2).Reduce(valkey2,cub::Max());
			//reduce on maximum count of shared hash table
			bool termination = false;
			if (threadIdx.x == 0) {
				//get the least 32-bit
				int max_sh_key = aggregate2;
				int max_sh_count = aggregate2>>32;
				int max_cm_count = aggregate>>32;
				// if shared count > the maximum in CM sketch.
				if(max_sh_count > max_cm_count){
					termination = true;

					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						temp_labels[v] = my_max_key - 1;

					}
				}
			}

			bool term = __shfl_sync(FULL_MASK,termination,0);
			if(term == false){
				// if cm is invaild, read all labels again
				for (int i = begin; i < end; i += NT) {
					int j = i + threadIdx.x;
					int key = 0;
					if (j < end) {
						key = adj_labels[j]+1;
					}
					auto lmask = __match_any_sync(FULL_MASK,key);
					auto lleader  = __ffs(lmask) - 1;
					auto count  = __popc(lmask);

					if (g.thread_rank()== lleader&&key!=0){
						// check shared hash table, if not exists write to global hash table.
						auto hashkey = hash(key);
						hashkey = (hashkey)%W;
						auto old = keys[hashkey];
						if(old != key){

							atomicAdd(counter, 1);
							count = gtable.increment(begin,end,key,count);
							if(count >my_max_count){
								my_max_count = count;
								my_max_key   = key;
							}
						}
					}
				}
				//block reduce
				__syncthreads();
				typedef cub::BlockReduce<unsigned long long,NT> BlockReduce3;
				__shared__ typename BlockReduce3::TempStorage temp_storage3;
				unsigned long long valkey3 = ((unsigned long long)my_max_count <<32 )+ my_max_key;
				unsigned long long  aggregate3 = BlockReduce3(temp_storage3).Reduce(valkey3,cub::Max());
				if (threadIdx.x == 0) {
					//get the least 32-bit
					int my_max_key = aggregate3;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						temp_labels[v] = my_max_key - 1;
					}
				}
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn2
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable)
{
	__shared__  int vals[W];
    __shared__  int keys[W];
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;
				if (g.thread_rank()== lleader&&key!=0){
					//first hash to CM then check shared hash table
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0 || old == key){
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							int tempcount = atomicAdd(&vals[hashkey],count);
							count +=tempcount;
							if(my_max_count < count){
								my_max_count = count;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){
						count = gtable.increment(begin,end,key,count);

						if(count >my_max_count){
							my_max_count = count;
							my_max_key   = key;
						}
					}

				}
			}
			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
//					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}


template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn2_labelload
(int offset, int n, int *neighbors,int *offsets, int *labels_write,int *labels_read,int* vertices,int *counter, GlobalHT gtable)
{
	__shared__  int vals[W];
    __shared__  int keys[W];
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = labels_read[neighbors[j]]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;
				if (g.thread_rank()== lleader&&key!=0){
					//first hash to CM then check shared hash table
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0 || old == key){
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							int tempcount = atomicAdd(&vals[hashkey],count);
							count +=tempcount;
							if(my_max_count < count){
								my_max_count = count;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){
						count = gtable.increment(begin,end,key,count);

						if(count >my_max_count){
							my_max_count = count;
							my_max_key   = key;
						}
					}

				}
			}
			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = labels_read[v];
				if (lbl != my_max_key - 1) {
//					atomicAdd(counter, 1);
					labels_write[v] = my_max_key - 1;
				}
			}
      }
   }
}


template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn2_seed
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable, bool* seednodes)
{
	__shared__  int vals[W];
    __shared__  int keys[W];
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;

    		bool seed = seednodes[v];
			int my_max_key = v+1;
			int my_max_count = 0;
			if(seed == false){
				const int begin = offsets[v];
				const int end = offsets[v+1];

				//if there is single node, key = my_max_key.

				int inserted = 0;
				for (int i = begin; i < end; i += NT) {
					int j = i + threadIdx.x;
					int key = 0;
					if (j < end) {
						key = adj_labels[j]+1;
					}
					auto lmask = __match_any_sync(FULL_MASK,key);
					auto lleader  = __ffs(lmask) - 1;
					auto count  = __popc(lmask);
					bool hashToGlobalMemory = true;
					if (g.thread_rank()== lleader&&key!=0){
						//first hash to CM then check shared hash table
						auto hashkey = hash(key);
						hashkey = hashkey%W;
						auto old = keys[hashkey];
						if(old == 0 || old == key){
							old = atomicCAS(&keys[hashkey],0,key);
							if(old ==0|| old == key){
								int tempcount = atomicAdd(&vals[hashkey],count);
								count +=tempcount;
								if(my_max_count < count){
									my_max_count = count;
									my_max_key = key;
								}
								hashToGlobalMemory = false;
							}
						}
						if(hashToGlobalMemory == true){
							count = gtable.increment(begin,end,key,count);

							if(count >my_max_count){
								my_max_count = count;
								my_max_key   = key;
							}
						}

					}
				}
				//block reduce
				__syncthreads();
				typedef cub::BlockReduce< long long,NT> BlockReduce;
				__shared__ typename BlockReduce::TempStorage temp_storage;
				unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
				unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());
				if (threadIdx.x == 0) {
					//get the least 32-bit
					int my_max_key = aggregate;

					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
	//					atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
					}
				}
			}else{
				temp_labels[v] = my_max_key - 1;
			}

      }
   }
}


template<int TW, int VT,int NT>
__global__ void l_big_update_syn_no_shared
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable)
{

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){

//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				if (g.thread_rank()== lleader&&key!=0){

					count = gtable.increment(begin,end,key,count);

					if(count >my_max_count){
						my_max_count = count;
						my_max_key   = key;
					}
				}
			}
			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}
template<int TW, int VT,int NT>
__global__ void sl_big_update_syn_no_shared
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable,int* gcounter_w)
{

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){

//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				if (g.thread_rank()== lleader&&key!=0){

					count = gtable.increment(begin,end,key,count);

					if(count >my_max_count){
						my_max_count = count;
						my_max_key   = key;
					}
				}
			}
			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];
				gcounter_w[lbl-1] ++;
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}
template<int TW, int VT,int NT>
__global__ void ml_big_update_syn_no_shared
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable,int* gcounter_r, int* gcounter_w,int delta)
{

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){

//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_score = INT_MIN;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				int score = INT_MIN;
				if (g.thread_rank()== lleader&&key!=0){

					count = gtable.increment(begin,end,key,count);
					score = delta*count - gcounter_r[key-1];
					if(score >my_max_score){
						my_max_score = score;
						my_max_key   = key;
					}
				}
			}
//			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_score <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0 && v < n) {
				//get the least 32-bit
				int my_max_key = aggregate;
				auto lbl = temp_labels[v];

				if(lbl -1 >= 0){
					gcounter_w[lbl-1]++;
				}

				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn4
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable)
{
	__shared__  int vals[W];
    __shared__  int keys[W];
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;
				if (g.thread_rank()== lleader&&key!=0){
					//first hash to CM then check shared hash table
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0 || old == key){
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							int tempcount = atomicAdd(&vals[hashkey],count);
							count +=tempcount;
							if(my_max_count < count){
								my_max_count = count;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
				}
				auto ghashkey = hash(key)%(end-begin) + begin;
				auto hmask = __ballot_sync(FULL_MASK, hashToGlobalMemory);
				while(hmask != 0){
					auto leader = __ffs(hmask)-1;
					hmask &= !(1<<leader);
					for(int p = begin;p<end;p+=32){
						auto q = p + g.thread_rank();
						// mark check_key = -1 as invaild value
						uint32_t check_key = -1;
						if(q<end){
							check_key = gtable.keys[q];
						}
						//check for key
						auto maskforkey = __ballot_sync(FULL_MASK,key==check_key);
						auto hashcount = 0;
						//if find
						if(maskforkey!=0){
							if(g.thread_rank() == __ffs(maskforkey)-1){
								hashcount = atomicAdd(&(gtable.vals[q]),count);
								hashcount += count;
							}
							if(hashcount > my_max_count){
								my_max_count = hashcount;
								my_max_key = key;
							}
							break;
						}
						//check for zero and insert
						maskforkey = __ballot_sync(FULL_MASK,check_key==0);
						if(maskforkey!=0){
							if(g.thread_rank() == __ffs(maskforkey)-1)
							{
								int old = atomicCAS(&(gtable.keys[q]),0,key);
								if(old == 0){
									hashcount = atomicAdd(&(gtable.vals[q]),count);
									hashcount += count;
									if(hashcount > my_max_count){
										my_max_count = hashcount;
										my_max_key = key;
									}
									break;
								}
							}
						}
					}
				}
			}
			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn2
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable,bool* if_update)
{
	__shared__  int vals[W];
    __shared__  int keys[W];
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;
				if (g.thread_rank()== lleader&&key!=0){
					//first hash to CM then check shared hash table
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0 || old == key){
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							int tempcount = atomicAdd(&vals[hashkey],count);
							count +=tempcount;
							if(my_max_count < count){
								my_max_count = count;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){
						count = gtable.increment(begin,end,key,count);
						atomicAdd(counter, 1);
						if(count >my_max_count){
							my_max_count = count;
							my_max_key   = key;
						}
					}
				}
			}
			//block reduce
			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					if_update[v] = 1;
//					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void l_big_update_syn2_seed
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable,bool* if_update, bool* seednodes)
{
	__shared__  int vals[W];
    __shared__  int keys[W];
    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
    		bool seed = seednodes[v];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			if(!seed){
				int inserted = 0;
				const int begin = offsets[v];
				const int end = offsets[v+1];
				for (int i = begin; i < end; i += NT) {
					int j = i + threadIdx.x;
					int key = 0;
					if (j < end) {
						key = adj_labels[j]+1;
					}
					auto lmask = __match_any_sync(FULL_MASK,key);
					auto lleader  = __ffs(lmask) - 1;
					auto count  = __popc(lmask);
					bool hashToGlobalMemory = true;
					if (g.thread_rank()== lleader&&key!=0){
						//first hash to CM then check shared hash table
						auto hashkey = hash(key);
						hashkey = hashkey%W;
						auto old = keys[hashkey];
						if(old == 0 || old == key){
							old = atomicCAS(&keys[hashkey],0,key);
							if(old ==0|| old == key){
								int tempcount = atomicAdd(&vals[hashkey],count);
								count +=tempcount;
								if(my_max_count < count){
									my_max_count = count;
									my_max_key = key;
								}
								hashToGlobalMemory = false;
							}
						}
						if(hashToGlobalMemory == true){
							count = gtable.increment(begin,end,key,count);
							atomicAdd(counter, 1);
							if(count >my_max_count){
								my_max_count = count;
								my_max_key   = key;
							}
						}
					}
				}
				//block reduce
				__syncthreads();
				typedef cub::BlockReduce< long long,NT> BlockReduce;
				__shared__ typename BlockReduce::TempStorage temp_storage;
				unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
				unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

				if (threadIdx.x == 0) {
					//get the least 32-bit
					int my_max_key = aggregate;

					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						if_update[v] = 1;
	//					atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
					}
				}
			}else{
				if (threadIdx.x == 0) {
					temp_labels[v] = my_max_key - 1;
				}
			}
      }
   }
}


// W > distinct label number.
template<int TW, int VT,int NT,int W>
__global__ void ml_big_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable, int* gcounter,int* gcounter_temp, int delta)
{
	 __shared__  int vals[W];
     __shared__  int keys[W];

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}

    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = INT_MIN;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = adj_labels[j]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;

				if (g.thread_rank()== lleader&&key!=0){
					// first hash to shared memory, if collision, hash to global memory
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0){
						// empty location
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							atomicAdd(&vals[hashkey],count);
							int score = count*delta - gcounter[key-1];
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}else{
						if(keys[hashkey]==key){
							count += vals[hashkey];
							vals[hashkey] = count;
							int score = count *delta -gcounter[key-1];
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){

						count = gtable.increment(begin,end,key,count);
						int score = count *delta -gcounter[key-1];
						if(score >my_max_count){
							my_max_count = score;
							my_max_key   = key;
						}
					}
				}
			}

			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			 long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			 long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];

				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void ml_big_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable, int* gcounter,int* gcounter_temp, int delta, bool* if_update)
{
   __shared__  int vals[W];
   __shared__  int keys[W];

   auto block = cg::this_thread_block();
   auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
   for (int stride = 0 ;stride<VT;++stride){
   	if(VT*blockIdx.x + stride+ offset < n){
   		// clear shared hash table and global table
   		for(int stride = 0;stride<W ;stride+=NT ){
   			if(threadIdx.x+stride <W){
   				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
   			}
   		}

   		__syncthreads();
//		const int v = vertices[VT*blockIdx.x + stride+ offset];
   		const int v = VT*blockIdx.x + stride+ offset;
		const int begin = offsets[v];
		const int end = offsets[v+1];
		//if there is single node, key = my_max_key.
		int my_max_key = v+1;
		int my_max_count = INT_MIN;
		int inserted = 0;
		for (int i = begin; i < end; i += NT) {
			int j = i + threadIdx.x;
			int key;
			if (j < end) {
				key = adj_labels[j]+1;
			}else{
				key = 0;
			}
			auto lmask = __match_any_sync(FULL_MASK,key);
			auto lleader  = __ffs(lmask) - 1;
			auto count  = __popc(lmask);
			bool hashToGlobalMemory = true;

			if (g.thread_rank()== lleader&&key!=0){
				// first hash to shared memory, if collision, hash to global memory
				auto hashkey = hash(key);
				hashkey = hashkey%W;
				auto old = keys[hashkey];

				if(old == 0){
					// empty location
					old = atomicCAS(&keys[hashkey],0,key);
					if(old ==0|| old == key){
						atomicAdd(&vals[hashkey],count);
						int score = count*delta - gcounter[key-1];
						if(my_max_count < score){
							my_max_count = score;
							my_max_key = key;
						}
						hashToGlobalMemory = false;
					}
				}else{
					if(keys[hashkey]==key){
						count += vals[hashkey];
						vals[hashkey] = count;
						int score = count *delta -gcounter[key-1];
						if(my_max_count < score){
							my_max_count = score;
							my_max_key = key;
						}
						hashToGlobalMemory = false;
					}
				}
				if(hashToGlobalMemory == true){

					count = gtable.increment(begin,end,key,count);
					int score = count *delta -gcounter[key-1];
					if(score >my_max_count){
						my_max_count = score;
						my_max_key   = key;
					}
				}
			}
		}

		__syncthreads();
		typedef cub::BlockReduce< long long,NT> BlockReduce;
		__shared__ typename BlockReduce::TempStorage temp_storage;
		long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
		long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

		if (threadIdx.x == 0) {
			//get the least 32-bit
			int my_max_key = aggregate;

			auto lbl = temp_labels[v];

			if (lbl != my_max_key - 1) {
				atomicAdd(counter, 1);
				temp_labels[v] = my_max_key - 1;
				if_update[v] = 1;
			}
			atomicAdd(&gcounter_temp[my_max_key - 1],1);
		}
     }
  }
}

template<int TW, int VT,int NT,int W>
__global__ void ml_big_update_syn_labelload
(int offset, int n, int *neighbors,int *offsets, int *labels_write,int *labels_read,int* vertices,int *counter, GlobalHT gtable, int* gcounter,int* gcounter_temp, int delta)
{
	 __shared__  int vals[W];
     __shared__  int keys[W];

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}

    		__syncthreads();
//			const int v = vertices[VT*blockIdx.x + stride+ offset];
    		const int v = VT*blockIdx.x + stride+ offset;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = INT_MIN;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key = 0;
				if (j < end) {
					key = labels_read[neighbors[j]]+1;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;

				if (g.thread_rank()== lleader&&key!=0){
					// first hash to shared memory, if collision, hash to global memory
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0){
						// empty location
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							atomicAdd(&vals[hashkey],count);
							int score = count*delta - gcounter[key-1];
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}else{
						if(keys[hashkey]==key){
							count += vals[hashkey];
							vals[hashkey] = count;
							int score = count *delta -gcounter[key-1];
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){

						count = gtable.increment(begin,end,key,count);
						int score = count *delta -gcounter[key-1];
						if(score >my_max_count){
							my_max_count = score;
							my_max_key   = key;
						}
					}
				}
			}

			__syncthreads();
			typedef cub::BlockReduce< long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			 long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			 long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = labels_read[v];

				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					labels_write[v] = my_max_key - 1;
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void ml_big_update_syn_opt
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable, int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha)
{
	 __shared__  int vals[W];
     __shared__  int keys[W];

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}

    		__syncthreads();
			const int v = vertices[VT*blockIdx.x + stride+ offset];
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = INT_MIN;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key;
				if (j < end) {
					key = adj_labels[j]+1;
				}else{
					key = 0;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;
				if (g.thread_rank()== lleader&&key!=0){
					// first hash to shared memory, if collision, hash to global memory
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0){
						// empty location
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							atomicAdd(&vals[hashkey],count);
							int score = count*delta - (gcounterlist[j]+alpha);
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}else{
						if(keys[hashkey]==key){
							count += vals[hashkey];
							vals[hashkey] = count;
							int score = count*delta - (gcounterlist[j]+alpha);
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){

						count = gtable.increment(begin,end,key,count);
						int score = count*delta - (gcounterlist[j]+alpha);
						if(score >my_max_count){
							my_max_count = score;
							my_max_key   = key;
						}
					}
				}
			}
			// find lower bound.
			typedef cub::BlockReduce<int,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			int lower_count = BlockReduce(temp_storage).Reduce(my_max_count,cub::Max());
			int lowercount = __shfl_sync(FULL_MASK,lower_count,0);
			//
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key;
				if (j < end) {
					key = adj_labels[j]+1;
				}else{
					key = 0;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				bool hashToGlobalMemory = true;
				int count = 0;
				if (g.thread_rank()== lleader&&key!=0){
					// first hash to shared memory, if collision, hash to global memory
					auto hashkey = hash(key);
					hashkey = hashkey%W;

					if(keys[hashkey]==key){
						count = vals[hashkey];
						hashToGlobalMemory = false;
					}
					if(hashToGlobalMemory == true){
						count = gtable.search(begin,end,key);
					}
					int score = count*delta - (gcounterlist[j]-alpha);
					if(score >lowercount){
						score = count*delta - gcounter[j];
						if(score >my_max_count){
							my_max_count = score;
							my_max_key   = key;
						}
					}
				}
			}


			typedef cub::BlockReduce< long long,NT> BlockReduce2;
			__shared__ typename BlockReduce2::TempStorage temp_storage2;
			 long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			 long long aggregate = BlockReduce2(temp_storage2).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];

				if (lbl != my_max_key - 1) {

					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}
      }
   }
}

template<int TW, int VT,int NT,int W>
__global__ void ml_big_update_syn_opt
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, GlobalHT gtable, int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha, bool* if_update)
{
	 __shared__  int vals[W];
     __shared__  int keys[W];

    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		// clear shared hash table and global table
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){
    				vals[threadIdx.x+stride] = 0;
					keys[threadIdx.x+stride] = 0;
    			}
    		}

    		__syncthreads();
			const int v = vertices[VT*blockIdx.x + stride+ offset];
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = INT_MIN;
			int inserted = 0;
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key;
				if (j < end) {
					key = adj_labels[j]+1;
				}else{
					key = 0;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);
				bool hashToGlobalMemory = true;
				if (g.thread_rank()== lleader&&key!=0){
					// first hash to shared memory, if collision, hash to global memory
					auto hashkey = hash(key);
					hashkey = hashkey%W;
					auto old = keys[hashkey];
					if(old == 0){
						// empty location
						old = atomicCAS(&keys[hashkey],0,key);
						if(old ==0|| old == key){
							atomicAdd(&vals[hashkey],count);
							int score = count*delta - (gcounterlist[j]+alpha);
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}

					}else{
						if(keys[hashkey]==key){
							count += vals[hashkey];
							vals[hashkey] = count;
							int score = count*delta - (gcounterlist[j]+alpha);
							if(my_max_count < score){
								my_max_count = score;
								my_max_key = key;
							}
							hashToGlobalMemory = false;
						}
					}
					if(hashToGlobalMemory == true){

						count = gtable.increment(begin,end,key,count);

						int score = count*delta - (gcounterlist[j]+alpha);
						if(score >my_max_count){
							my_max_count = score;
							my_max_key   = key;
						}
					}
				}
			}
//			// find lower bound.
			typedef cub::BlockReduce<int,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			int lower_count = BlockReduce(temp_storage).Reduce(my_max_count,cub::Max());
			int lowercount = __shfl_sync(FULL_MASK,lower_count,0);
//
			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key;
				if (j < end) {
					key = adj_labels[j]+1;
				}else{
					key = 0;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				bool hashToGlobalMemory = true;
				int count = 0;
				if (g.thread_rank()== lleader&&key!=0){
					// first hash to shared memory, if collision, hash to global memory
					auto hashkey = hash(key);
					hashkey = hashkey%W;

					if(keys[hashkey]==key){
						count = vals[hashkey];
						hashToGlobalMemory = false;
					}
					if(hashToGlobalMemory == true){
						count = gtable.search(begin,end,key);
					}

					int score = count*delta - (gcounterlist[j]-alpha);
					if(score >lowercount){
						score = count*delta - gcounter[key-1];
						if(score >my_max_count){
							my_max_count = score;
							my_max_key   = key;
						}
					}
				}
			}
			__syncthreads();
//
			typedef cub::BlockReduce< long long,NT> BlockReduce2;
			__shared__ typename BlockReduce2::TempStorage temp_storage2;
			 long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			 long long aggregate = BlockReduce2(temp_storage2).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];

				if (lbl != my_max_key - 1) {
					if_update[v] = 1;
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}
      }
   }
}



// my method one warp one node.
template<int TW, int VT,int NT,int BW,int W>
__global__ void l_big_update_syn_p
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,double *prob,double *rand_prob)
{

	 __shared__ int h1vals[W];
	 __shared__ int h2vals[W];
	 __shared__ int h3vals[W];


    auto block = cg::this_thread_block();
    auto g = cg::tiled_partition<TW>(block);

	#pragma unroll
    for (int stride = 0 ;stride<VT;++stride){
    	if(VT*blockIdx.x + stride+ offset < n){
    		for(int stride = 0;stride<W ;stride+=NT ){
    			if(threadIdx.x+stride <W){

    				h1vals[threadIdx.x+stride] = 0;
					h2vals[threadIdx.x+stride] = 0;
					h3vals[threadIdx.x+stride] = 0;
    			}
    		}
    		__syncthreads();
			const int v = vertices[VT*blockIdx.x + stride+ offset];
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//if there is single node, key = my_max_key.
			int my_max_key = v+1;
			int my_max_count = 0;
			int inserted = 0;
			double p = 1-prob[v];
			if(rand_prob[v]>p){
				continue;
			}

			for (int i = begin; i < end; i += NT) {
				int j = i + threadIdx.x;
				int key;
				if (j < end) {
					key = adj_labels[j]+1;
				}else{
					key = 0;
				}
				auto lmask = __match_any_sync(FULL_MASK,key);
				auto lleader  = __ffs(lmask) - 1;
				auto count  = __popc(lmask);

				if (g.thread_rank()== lleader&&key!=0){
					int h1 = hash1(key)%W;
					int h2 = hash2(key)%W;
					int h3 = hash3(key)%W;
					int mincount;
					int count1 = atomicAdd(&h1vals[h1],count)+count;
					int count2 = atomicAdd(&h2vals[h2],count)+count;
					int count3 = atomicAdd(&h3vals[h3],count)+count;
					mincount = min(count1,count2);
					mincount = min(mincount,count3);

					if(mincount >my_max_count){
						my_max_count = mincount;
						my_max_key   = key;
					}
				}
			}

			__syncthreads();
			typedef cub::BlockReduce<unsigned long long,NT> BlockReduce;
			__shared__ typename BlockReduce::TempStorage temp_storage;
			unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			unsigned long long  aggregate = BlockReduce(temp_storage).Reduce(valkey,cub::Max());

			if (threadIdx.x == 0) {
				//get the least 32-bit
				int my_max_key = aggregate;

				auto lbl = temp_labels[v];

				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
					prob[v] = 0;
				}else{
					prob[v] = 0.5+prob[v]/2;
				}
			}
      }
   }
}





// my method one warp one node.
template<int TW,int VT, int NT, int d>
__global__ void l_medium_update_syn
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);

	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = 0;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
//			int v = vertices[id];
			const int v = id;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }
//				 if(v == 443511){
//					 printf("lane id %d,j %d, mylanekey %d,v %d,neighbors %d \n",g.thread_rank(),j,mylanekey,v,neighbors[j]);
//				 }
				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;

								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);


			 typedef cub::WarpReduce<unsigned long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
//					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
			}

		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void l_medium_update_syn_labelload
(int offset,int n, int *neighbors,int *offsets, int *labels_write,int *labels_read,int* vertices,int *counter,const int medium)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);

	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = 0;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
//			int v = vertices[id];
			const int v = id;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = labels_read[neighbors[j]]+1;
				 }else{
					 mylanekey = 0;
				 }

				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;

								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);


			 typedef cub::WarpReduce<unsigned long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = labels_read[v];
				if (lbl != my_max_key - 1) {
//					atomicAdd(counter, 1);
					labels_write[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
			}

		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void l_medium_update_syn_seed
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium, bool* seednodes)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);

	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = 0;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
//			int v = vertices[id];
			const int v = id;
			bool seed = seednodes[v];
			if(!seed){
				const int begin = offsets[v];
				const int end = offsets[v+1];
				//load value
				// number of neighbors is medium;
				for (int i = begin; i < end; i += TW) {
					int j = i + g.thread_rank();
					 // sample first TW label
					 if(j< end){
						 mylanekey = adj_labels[j]+1;
					 }else{
						 mylanekey = 0;
					 }
	//				 if(v == 443511){
	//					 printf("lane id %d,j %d, mylanekey %d,v %d,neighbors %d \n",g.thread_rank(),j,mylanekey,v,neighbors[j]);
	//				 }
					 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
					 auto lleader  = __ffs(lmask) - 1;
					 auto count = __popc(lmask);
					 if( mylanekey!=0 && lleader==g.thread_rank()){
						auto hashkey = hash(mylanekey);
						#pragma unroll
						for(int k = 0;k<medium;++k){
							hashkey = (hashkey+k)%medium;
							auto old = keys[hashkey + warpid*medium];
							if(old == 0){
								// empty location
								old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
								if(old ==0){
									vals[hashkey + warpid*medium] = count;
									if(my_max_count < count){
										my_max_count = count;
										my_max_key = mylanekey;
									}
									break;
								}
							}else{
								if(keys[hashkey + warpid*medium]==mylanekey){
									count += vals[hashkey + warpid*medium];
									vals[hashkey + warpid*medium] = count;

									if(my_max_count < count){
										my_max_count = count;
										my_max_key = mylanekey;
									}
									break;
								}
							}
						}

					 }
				 }
				 __syncwarp(FULL_MASK);


				 typedef cub::WarpReduce<unsigned long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
				 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());

				if(g.thread_rank()==0){
					my_max_key = aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
	//					atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}else{
				if(g.thread_rank()==0){
					temp_labels[v] = v;
				}
			}
		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void l_medium_update_syn
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium, bool* if_update)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = 0;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
//			int v = vertices[id];
			const int v= id;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;

								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);
			 typedef cub::WarpReduce<unsigned long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
//					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
					if_update[v] = 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
			}

		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void l_medium_update_syn_seed
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium, bool* if_update, bool* seednodes)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = 0;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
//			int v = vertices[id];
			const int v= id;
			bool seed = seednodes[v];
			if(!seed){
				const int begin = offsets[v];
				const int end = offsets[v+1];
				//load value
				// number of neighbors is medium;
				for (int i = begin; i < end; i += TW) {
					int j = i + g.thread_rank();
					 // sample first TW label
					 if(j< end){
						 mylanekey = adj_labels[j]+1;
					 }else{
						 mylanekey = 0;
					 }

					 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
					 auto lleader  = __ffs(lmask) - 1;
					 auto count = __popc(lmask);
					 if( mylanekey!=0 && lleader==g.thread_rank()){
						auto hashkey = hash(mylanekey);
						#pragma unroll
						for(int k = 0;k<medium;++k){
							hashkey = (hashkey+k)%medium;
							auto old = keys[hashkey + warpid*medium];
							if(old == 0){
								// empty location
								old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
								if(old ==0){
									vals[hashkey + warpid*medium] = count;
									if(my_max_count < count){
										my_max_count = count;
										my_max_key = mylanekey;
									}
									break;
								}
							}else{
								if(keys[hashkey + warpid*medium]==mylanekey){
									count += vals[hashkey + warpid*medium];
									vals[hashkey + warpid*medium] = count;

									if(my_max_count < count){
										my_max_count = count;
										my_max_key = mylanekey;
									}
									break;
								}
							}
						}

					 }
				 }
				 __syncwarp(FULL_MASK);
				 typedef cub::WarpReduce<unsigned long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
				 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());

				if(g.thread_rank()==0){
					my_max_key = aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
	//					atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
						if_update[v] = 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}else{
				temp_labels[v] = v;
			}
		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void ml_medium_update_syn
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium, int* gcounter,int* gcounter_temp, int delta)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = INT_MIN;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }
//				 if(v == 443511){
//					 printf("lane id %d,j %d, mylanekey %d,v %d,neighbors %d \n",g.thread_rank(),j,mylanekey,v,neighbors[j]);
//				 }
				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - gcounter[mylanekey-1];
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - gcounter[mylanekey-1];
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);
			 typedef cub::WarpReduce< long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			  long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			  long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);

					temp_labels[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
//				if(my_max_key - 1< 0){
//					printf("!!!!!!!!!!!!!! \n");
//				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}

		}
	}
}


template<int TW,int VT, int NT, int d>
__global__ void ml_medium_update_syn_labelload
(int offset,int n, int *neighbors,int *offsets, int *labels_write,int *labels_read,int* vertices,int *counter,const int medium, int* gcounter,int* gcounter_temp, int delta)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = INT_MIN;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = labels_read[neighbors[j]]+1;
				 }else{
					 mylanekey = 0;
				 }
//				 if(v == 443511){
//					 printf("lane id %d,j %d, mylanekey %d,v %d,neighbors %d \n",g.thread_rank(),j,mylanekey,v,neighbors[j]);
//				 }
				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - gcounter[mylanekey-1];
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - gcounter[mylanekey-1];
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);
			 typedef cub::WarpReduce< long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			  long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			  long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = labels_read[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);

					labels_write[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}

				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}

		}
	}
}
template<int TW,int VT, int NT, int d>
__global__ void ml_medium_update_syn
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium, int* gcounter,int* gcounter_temp, int delta, bool* if_update)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = INT_MIN;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }
//				 if(v == 443511){
//					 printf("lane id %d,j %d, mylanekey %d,v %d,neighbors %d \n",g.thread_rank(),j,mylanekey,v,neighbors[j]);
//				 }
				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - gcounter[mylanekey-1];
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - gcounter[mylanekey-1];
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);
			 typedef cub::WarpReduce< long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			  long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			  long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					if_update[v] = 1;
					temp_labels[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}

		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void ml_medium_update_syn_opt
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium,  int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = INT_MIN;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = vertices[id];
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - (gcounterlist[j]+alpha);
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - (gcounterlist[j]+alpha);
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }

			// use lower bound.
			typedef cub::WarpReduce<int> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			int lower_count = WarpReduce(temp_storage[warpid]).Reduce(
					 my_max_count, cub::Max());
			int lowercount = __shfl_sync(FULL_MASK,lower_count,0);

			//
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;

				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						if(keys[hashkey + warpid*medium]==mylanekey){
							int count = vals[hashkey + warpid*medium];

							int score = delta*count - (gcounterlist[j]-alpha);
							if(lowercount < score){
								score = delta * count - gcounter[mylanekey-1];
								my_max_count = score;
								my_max_key = mylanekey;
							}
							break;
						}

					}

				 }
			 }
			//
			 typedef cub::WarpReduce< long long> WarpReduce2;
			 __shared__ typename WarpReduce2::TempStorage temp_storage2[NT/TW];
			 long long valkey = (( long long)my_max_count <<32 )+ my_max_key;
			 long long aggregate = WarpReduce2(temp_storage2[warpid]).Reduce(valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}

		}
	}
}

template<int TW,int VT, int NT, int d>
__global__ void ml_medium_update_syn_opt
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium,  int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha,bool* if_update)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = INT_MIN;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = vertices[id];
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - (gcounterlist[j]+alpha);
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;
								int score = delta*count - (gcounterlist[j]+alpha);
								if(my_max_count < score){
									my_max_count = score;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }

			// use lower bound.
			typedef cub::WarpReduce<int > WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];

			int lower_count = WarpReduce(temp_storage[warpid]).Reduce(
					 my_max_count, cub::Max());
			int lowercount = __shfl_sync(FULL_MASK,lower_count,0);

			//
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;

				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						if(keys[hashkey + warpid*medium]==mylanekey){
							int count = vals[hashkey + warpid*medium];

							int score = delta*count - (gcounterlist[j]-alpha);
							if(lowercount < score){
								score = delta * count - gcounter[mylanekey-1];
								my_max_count = score;
								my_max_key = mylanekey;
							}
							break;
						}

					}

				 }
			 }
			//
			 typedef cub::WarpReduce<long long> WarpReduce2;
			 __shared__ typename WarpReduce2::TempStorage temp_storage2[NT/TW];
			 long long valkey = ((long long)my_max_count <<32 )+ my_max_key;
			 long long aggregate = WarpReduce2(temp_storage2[warpid]).Reduce(valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					if_update[v] = 1;
					temp_labels[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
				}
				atomicAdd(&gcounter_temp[my_max_key - 1],1);
			}

		}
	}
}

// my method one warp one node.
template<int TW,int VT, int NT, int d>
__global__ void l_medium_update_syn_p
(int offset,int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,const int medium,double *prob,double *rand_prob)
{
	 __shared__  int vals[NT*d];
	 __shared__  int keys[NT*d];
	int mylanekey;
	int my_max_key;
	int my_max_count;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		for (int i = 0;i< d;++i){
			vals[threadIdx.x+i*NT] =0;
			keys[threadIdx.x+i*NT] =0;
		}
		my_max_key = 0;
		my_max_count = 0;
		__syncthreads();
	//	 the index of warp(group)
		int warpid =  threadIdx.x / TW;
	//	 the vertex about to process
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = vertices[id];
			double p = 1-prob[v];
			if(rand_prob[v]>p){
				return;
			}
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			// number of neighbors is medium;
			for (int i = begin; i < end; i += TW) {
				int j = i + g.thread_rank();
				 // sample first TW label
				 if(j< end){
					 mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }
//				 if(v == 443511){
//					 printf("lane id %d,j %d, mylanekey %d,v %d,neighbors %d \n",g.thread_rank(),j,mylanekey,v,neighbors[j]);
//				 }
				 auto lmask = __match_any_sync(FULL_MASK,mylanekey);
				 auto lleader  = __ffs(lmask) - 1;
				 auto count = __popc(lmask);
				 if( mylanekey!=0 && lleader==g.thread_rank()){
					auto hashkey = hash(mylanekey);
					#pragma unroll
					for(int k = 0;k<medium;++k){
						hashkey = (hashkey+k)%medium;
						auto old = keys[hashkey + warpid*medium];
						if(old == 0){
							// empty location
							old = atomicCAS(&keys[hashkey + warpid*medium],0,mylanekey);
							if(old ==0){
								vals[hashkey + warpid*medium] = count;
								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}else{
							if(keys[hashkey + warpid*medium]==mylanekey){
								count += vals[hashkey + warpid*medium];
								vals[hashkey + warpid*medium] = count;

								if(my_max_count < count){
									my_max_count = count;
									my_max_key = mylanekey;
								}
								break;
							}
						}
					}

				 }
			 }
			 __syncwarp(FULL_MASK);
			 typedef cub::WarpReduce<unsigned long long> WarpReduce;
			 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
			 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
			 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
					 valkey, cub::Max());

			if(g.thread_rank()==0){
				my_max_key = aggregate;
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
					//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					prob[v] = 0;
				}else{
					prob[v] = 0.5+prob[v]/2;
				}
			}

		}
	}
}

template<int TW, int NT>
__global__ void l_medium_update_syn2
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int* warp_v,int* warp_id,int* max_count,int *counter,GlobalHT gtable,int hashstart,int Maxwarp)
{
	//record max in lanes
	int my_max_count =0;
	int my_max_key = 0;
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	//get warp id among all threads
	const int warpIdAll = blockIdx.x*(NT/TW)+threadIdx.x/TW;
	int v = -1;
	int begin = 0;
	int end = 0;
	int warpIdBlock = 0;
	if(warpIdAll< Maxwarp){
		//compute which node to process according to warpIdAll
		v= warp_v[warpIdAll];
		// compute the index of warps given the same v.
		warpIdBlock = warp_id[warpIdAll];
		begin = offsets[v];
		end = offsets[v+1];
	}else{
		// last some empty warps in last block
	}

	//assign threads
	int j = begin + TW * warpIdBlock + g.thread_rank();

	int key;
	if (j < end) {
		key = adj_labels[j]+1;
	}else{
		key = 0;
	}



	auto lmask = __match_any_sync(FULL_MASK,key);
	auto lleader  = __ffs(lmask) - 1;
	auto count  = __popc(lmask);
	if (g.thread_rank()== lleader&&key!=0){
		int mincount = gtable.increment(begin-hashstart,end-hashstart,key,count);

		int c =mincount;
		if(c >my_max_count){
			my_max_count = c;
			my_max_key   = key;
		}
	}
	//warp reduce
	//	 the index of warp(group)
	 int warpid =  threadIdx.x / TW;
	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
	 unsigned long long valkey = ((unsigned long long)my_max_count <<32 )+ my_max_key;
	 unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
			 valkey, cub::Max());
	__syncthreads();
//	 Try to update the label of the vertex a warp handles
	if (g.thread_rank()== 0 && v!=-1) {
		my_max_count = valkey>>32;
		my_max_key = valkey;
		// across block reduce
		auto ret = atomicMax(&max_count[v], my_max_count);
		if (ret < my_max_count ) {
			auto lbl = temp_labels[v];
			temp_labels[v] = my_max_key - 1;

			if (lbl != my_max_key - 1) {
				// May count redundantly
//					atomicAdd(counter, 1);
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void l_small_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
				}

				 unsigned long long valkey = ((unsigned long long)count <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce< unsigned long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				  unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
//						atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void l_small_update_syn_labelload
(int offset, int n, int *neighbors,int *offsets, int *labels_write,int *labels_read,int* vertices,int *counter)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				labels_write[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = labels_read[neighbors[j]]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
				}

				 unsigned long long valkey = ((unsigned long long)count <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce< unsigned long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				  unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = labels_read[v];
					if (lbl != my_max_key - 1) {
//						atomicAdd(counter, 1);
						labels_write[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void l_small_update_syn_seed
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,bool* seednodes)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			bool seed = seednodes[v];
			if(!seed ){
				//printf(" %d \n",id);
				const int begin = offsets[v];
				const int end = offsets[v+1];
				//load value
				//no neighbor
				if(begin==end){
					my_max_key = v+1;
					temp_labels[v] = my_max_key - 1;
				}else{
					// number of neighbors is small;
					int j = g.thread_rank()+begin;
					 // sample first TW label
					 if(j< end){
						mylanekey = adj_labels[j]+1;
					 }else{
						 mylanekey = 0;
					 }

					auto gmask = __match_any_sync(FULL_MASK,mylanekey);

					int count = 0;
					if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
						count = __popc(gmask);
					}

					 unsigned long long valkey = ((unsigned long long)count <<32 )+ mylanekey;

					 // cub reduce
					 typedef cub::WarpReduce< unsigned long long> WarpReduce;
					 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
					  unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
							 valkey, cub::Max());
					 if(g.thread_rank()==0){
						my_max_key =  aggregate;
						auto lbl = temp_labels[v];
						if (lbl != my_max_key - 1) {
	//						atomicAdd(counter, 1);
							temp_labels[v] = my_max_key - 1;
							//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
						}
					}
				}
			}else{
				temp_labels[v] = v;
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void l_small_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int* counter,bool* if_update)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = id;
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
				}

				 unsigned long long valkey = ((unsigned long long)count <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce< unsigned long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				  unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						if_update[v] = 1;
//						atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void l_small_update_syn_seed
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int* counter,bool* if_update,bool* seednodes)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = id;
			bool seed= seednodes[v];
			//printf(" %d \n",id);
			if(!seed){
				const int begin = offsets[v];
				const int end = offsets[v+1];
				//load value
				//no neighbor
				if(begin==end){
					my_max_key = v+1;
					temp_labels[v] = my_max_key - 1;
				}else{
					// number of neighbors is small;
					int j = g.thread_rank()+begin;
					 // sample first TW label
					 if(j< end){
						mylanekey = adj_labels[j]+1;
					 }else{
						 mylanekey = 0;
					 }

					auto gmask = __match_any_sync(FULL_MASK,mylanekey);

					int count = 0;
					if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
						count = __popc(gmask);
					}

					 unsigned long long valkey = ((unsigned long long)count <<32 )+ mylanekey;

					 // cub reduce
					 typedef cub::WarpReduce< unsigned long long> WarpReduce;
					 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
					  unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
							 valkey, cub::Max());
					 if(g.thread_rank()==0){
						my_max_key =  aggregate;
						auto lbl = temp_labels[v];
						if (lbl != my_max_key - 1) {
							if_update[v] = 1;
	//						atomicAdd(counter, 1);
							temp_labels[v] = my_max_key - 1;
							//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
						}
					}
				}
			}else{
				temp_labels[v] = v;
			}
		}
	}
}
//degree 32 -16
template<int TW,int VT, int NT>
__global__ void ml_small_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,int* gcounter,int* gcounter_temp, int delta)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				int score = INT_MIN;

				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
					score = count*delta - gcounter[mylanekey-1];
				}

				  long long valkey = (( long long)score <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce<  long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				   long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						atomicAdd(counter, 1);

						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}

					atomicAdd(&gcounter_temp[my_max_key - 1],1);
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void ml_small_update_syn_labelload
(int offset, int n, int *neighbors,int *offsets, int *labels_write,int *labels_read,int* vertices,int *counter,int* gcounter,int* gcounter_temp, int delta)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				labels_write[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = labels_read[neighbors[j]]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				int score = INT_MIN;

				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
					score = count*delta - gcounter[mylanekey-1];
				}

				  long long valkey = (( long long)score <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce<  long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				   long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = labels_read[v];
					if (lbl != my_max_key - 1) {
						atomicAdd(counter, 1);

						labels_write[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}

					atomicAdd(&gcounter_temp[my_max_key - 1],1);
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void ml_small_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,int* gcounter,int* gcounter_temp, int delta, bool* if_update)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			const int v = id;
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				int score = INT_MIN;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
					score = count*delta - gcounter[mylanekey-1];
				}

				long long valkey = (( long long)score <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce<long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				 long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						if_update[v] = 1;
						atomicAdd(counter, 1);

						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
					atomicAdd(&gcounter_temp[my_max_key - 1],1);
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void ml_small_update_syn_opt
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = vertices[id];
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				int score = INT_MIN;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
					score = count*delta - (gcounterlist[j] + alpha);
				}

				//lower bound
				 typedef cub::WarpReduce< int> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				  int lower_count = WarpReduce(temp_storage[warpid]).Reduce(
						 score, cub::Max());
				  int lowercount = __shfl_sync(FULL_MASK,lower_count,0);
				  score = count*delta - (gcounterlist[j] - alpha);
				  if(score > lowercount){
						  score = count*delta - gcounter[mylanekey-1];
				  }
				  long long valkey = (( long long)score <<32 )+ mylanekey;


				 // cub reduce
				 typedef cub::WarpReduce<  long long> WarpReduce2;
				 __shared__ typename WarpReduce2::TempStorage temp_storage2[NT/TW];
				   long long aggregate = WarpReduce2(temp_storage2[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						atomicAdd(counter, 1);
						atomicAdd(&gcounter_temp[my_max_key - 1],1);
						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void ml_small_update_syn_opt
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha,bool* if_update)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = vertices[id];
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				int score = INT_MIN;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
					score = count*delta - (gcounterlist[j] + alpha);
				}

				//lower bound
				 typedef cub::WarpReduce< int> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				  int lower_count = WarpReduce(temp_storage[warpid]).Reduce(
						 score, cub::Max());
				  int lowercount = __shfl_sync(FULL_MASK,lower_count,0);
				  if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					  score = count*delta - (gcounterlist[j] - alpha);
					  if(score > lowercount){
							  score = count*delta - gcounter[mylanekey-1];
					  }
				  }
				  long long valkey = (( long long)score <<32 )+ mylanekey;
				 // cub reduce
				  __syncthreads();
				 typedef cub::WarpReduce<  long long> WarpReduce2;
				 __shared__ typename WarpReduce2::TempStorage temp_storage2[NT/TW];
				   long long aggregate = WarpReduce2(temp_storage2[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						atomicAdd(counter, 1);
						atomicAdd(&gcounter_temp[my_max_key - 1],1);
						if_update[v] = 1;
						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
					}
				}
			}
		}
	}
}

template<int TW,int VT, int NT>
__global__ void l_small_update_syn_p
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,double *prob,double *rand_prob)
{
	int mylanekey;
	int my_max_key;

	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<TW>(block);
	// the index of warp(group)
	int warpid =  threadIdx.x / TW;
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)/TW+ offset;
		if(id< n){
			int v = vertices[id];
			double p = 1-prob[v];
			if(rand_prob[v]>p){
				return;
			}
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is small;
				int j = g.thread_rank()+begin;
				 // sample first TW label
				 if(j< end){
					mylanekey = adj_labels[j]+1;
				 }else{
					 mylanekey = 0;
				 }

				auto gmask = __match_any_sync(FULL_MASK,mylanekey);

				int count = 0;
				if(g.thread_rank()== __ffs(gmask)-1 && mylanekey!=0){
					count = __popc(gmask);
				}

				 unsigned long long valkey = ((unsigned long long)count <<32 )+ mylanekey;

				 // cub reduce
				 typedef cub::WarpReduce< unsigned long long> WarpReduce;
				 __shared__ typename WarpReduce::TempStorage temp_storage[NT/TW];
				  unsigned long long aggregate = WarpReduce(temp_storage[warpid]).Reduce(
						 valkey, cub::Max());
				 if(g.thread_rank()==0){
					my_max_key =  aggregate;
					auto lbl = temp_labels[v];
					if (lbl != my_max_key - 1) {
						atomicAdd(counter, 1);
						temp_labels[v] = my_max_key - 1;
						//printf("change v %d, label %d to %d \n",v ,lbl,my_max_key - 1);
						prob[v] = 0;
					}else{
						prob[v] = 0.5+prob[v]/2;
					}
				}
			}
		}
	}
}

template<int VT,int W>
__global__ void l_tiny_update_syn
(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter)
{
	int mylanekey;
	int my_max_key;
	int vals[W];
	int keys[W];
	#pragma unroll
	for(int k = 0;k < W;++k){
		keys[k] = 0;
		vals[k] = 0;
	}
	// the vertex about to process
	#pragma unroll
	for (int stride = 0;stride < VT;++stride){
		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)+ offset;
		if(id< n){
			int my_max_count = 0;
			int v = vertices[id];
			//printf(" %d \n",id);
			const int begin = offsets[v];
			const int end = offsets[v+1];
			//load value
			//no neighbor
			if(begin==end){
				my_max_key = v+1;
				temp_labels[v] = my_max_key - 1;
			}else{
				// number of neighbors is tiny;
				int count = 0;

				for(int j = begin;j<end;j++){
					mylanekey = adj_labels[j]+1;
					if(v== 3949564){
						printf("label1 %d \n",mylanekey);
					}
					#pragma unroll
					for(int k = 0;k < W;++k){
						if(keys[k]== 0){
							keys[k] = mylanekey;
							vals[k] = 1;
							count =1;
							if(v== 3949564)
							printf("insert: %d, count %d\n",mylanekey,vals[k]);
							break;
						}else{
							if(keys[k]==mylanekey){
								count = vals[k]+1;
//								vals[k] = count;
								vals[k]+=1;
								if(v== 3949564)
								printf("add: %d, count %d, vals[k] %d\n",mylanekey, count,vals[k]);
								break;
							}
						}
					}

					if(count > my_max_count){
						my_max_count = count;
						my_max_key = mylanekey;

					}
				}
				auto lbl = temp_labels[v];
				if (lbl != my_max_key - 1) {
					atomicAdd(counter, 1);
					temp_labels[v] = my_max_key - 1;
					if(v== 3949564){
						printf("Mcount1 %d \n",my_max_count);
					}
				}
			}
		}
	}
}

__global__ void l_tiny_update_syn2
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
	 }

	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 unsigned long long valkey = ((unsigned long long)count <<32 )+ label;
	 unsigned long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
//			atomicAdd(counter, 1);
			 temp_labels[v] = label - 1;
		 }
	 }
}

__global__ void l_tiny_update_syn2_labelload
( int *neighbors,int *offsets, int *labels_write,int *labels_read,int *counter,int* warp_v,int* warp_begin,int warpnumber)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = labels_read[neighbors[begin+g.thread_rank()]]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
	 }

	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 unsigned long long valkey = ((unsigned long long)count <<32 )+ label;
	 unsigned long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = labels_read[v];
		 if (lbl != label - 1) {
//			atomicAdd(counter, 1);
			labels_write[v] = label - 1;
		 }
	 }
}

__global__ void l_tiny_update_syn2_seed
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber, bool* seednodes)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
	 }

	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 unsigned long long valkey = ((unsigned long long)count <<32 )+ label;
	 unsigned long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 if(!seednodes[v]){
			 label =  aggregate;
			 auto lbl = temp_labels[v];
			 if (lbl != label - 1) {
	//			atomicAdd(counter, 1);
				temp_labels[v] = label - 1;
			 }
		 }else{
			 temp_labels[v] = v;
		 }

	 }
}

__global__ void l_tiny_update_syn2
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,bool* if_update)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
	 }

	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 unsigned long long valkey = ((unsigned long long)count <<32 )+ label;
	 unsigned long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
			 if_update[v] =1;
//			atomicAdd(counter, 1);
			temp_labels[v] = label - 1;
		 }
	 }
}

__global__ void l_tiny_update_syn2_seed
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,bool* if_update, bool* seednodes)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
	 }

	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 unsigned long long valkey = ((unsigned long long)count <<32 )+ label;
	 unsigned long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 if(!seednodes[v]){
			 label =  aggregate;
			 auto lbl = temp_labels[v];
			 if (lbl != label - 1) {
				 if_update[v] =1;
	//			atomicAdd(counter, 1);
				temp_labels[v] = label - 1;
			 }
		 }else{
			 temp_labels[v] = v;
		 }
	 }
}

__global__ void ml_tiny_update_syn
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,int* gcounter,int* gcounter_temp, int delta)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int score = INT_MIN;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }


	 if(g.thread_rank()== __ffs(lmask)-1&& v!=-1){
		 count = __popc(lmask);
		 score = delta*count - gcounter[label-1];
	 }

	 typedef cub::WarpReduce< long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 long long valkey = (( long long)score <<32 )+ label;
	 long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
			atomicAdd(counter, 1);
			temp_labels[v] = label - 1;
		 }
		 atomicAdd(&gcounter[label-1],1);
	 }
}

__global__ void ml_tiny_update_syn_labelload
( int *neighbors,int *offsets, int *labels_write,int *labels_read,int *counter,int* warp_v,int* warp_begin,int warpnumber,int* gcounter,int* gcounter_temp, int delta)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = labels_read[neighbors[begin+g.thread_rank()]]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int score = INT_MIN;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }


	 if(g.thread_rank()== __ffs(lmask)-1&& v!=-1){
		 count = __popc(lmask);
		 score = delta*count - gcounter[label-1];
	 }

	 typedef cub::WarpReduce< long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 long long valkey = (( long long)score <<32 )+ label;
	 long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = labels_read[v];
		 if (lbl != label - 1) {
			atomicAdd(counter, 1);
			labels_write[v] = label - 1;
		 }
		 atomicAdd(&gcounter[label-1],1);
	 }
}

__global__ void ml_tiny_update_syn
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,int* gcounter,int* gcounter_temp, int delta, bool * if_update)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int score = INT_MIN;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1&& v!=-1){
		 count = __popc(lmask);
		 score = delta*count - gcounter[label-1];
	 }

	 typedef cub::WarpReduce< long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	  long long valkey = (( long long)score <<32 )+ label;
	  long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
			atomicAdd(counter, 1);
			if_update[v] = 1;
			temp_labels[v] = label - 1;
		 }
		 atomicAdd(&gcounter[label-1],1);
	 }
}

__global__ void ml_tiny_update_syn_opt
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 int score = INT_MIN;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1&& v!=-1){
		 count = __popc(lmask);
		 score = count*delta - gcounterlist[begin+g.thread_rank()];
	 }
	 typedef cub::WarpReduce<int> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;

	 int lower_count = WarpReduce(temp_storage).HeadSegmentedReduce(
			 score, head, cub::Max());
	 int lowercount = __shfl_sync(vmask,lower_count,__ffs(vmask)-1);

	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
		 score = count*delta - (gcounterlist[begin-g.thread_rank()]+alpha);
		 if(score > lowercount ){
			 score = count*delta - (gcounterlist[begin-g.thread_rank()]-alpha);
		 }
	 }
	 typedef cub::WarpReduce< long long> WarpReduce2;
	 __shared__ typename WarpReduce2::TempStorage temp_storage2;
	 long long valkey = ((long long)count <<32 )+ label;
	 long long aggregate = WarpReduce2(temp_storage2).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
			atomicAdd(&gcounter[label-1],1);
			atomicAdd(counter, 1);
			temp_labels[v] = label - 1;
		 }
	 }
}

__global__ void ml_tiny_update_syn_opt
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,int* gcounter,int* gcounter_temp,int* gcounterlist, int delta, int alpha, bool* if_update)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];

	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 int score = INT_MIN;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
		 score = count*delta - gcounterlist[begin+g.thread_rank()];
	 }
	 typedef cub::WarpReduce<int> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;

	 int lower_count = WarpReduce(temp_storage).HeadSegmentedReduce(
			 score, head, cub::Max());
	 int lowercount = __shfl_sync(vmask,lower_count,__ffs(vmask)-1);

	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
		 score = count*delta - (gcounterlist[begin-g.thread_rank()]+alpha);
		 if(score > lowercount ){
			 score = count*delta - (gcounter[label-1]);
		 }
	 }
	 typedef cub::WarpReduce< long long> WarpReduce2;
	 __shared__ typename WarpReduce2::TempStorage temp_storage2;
	 long long valkey = ((long long)count <<32 )+ label;
	 long long aggregate = WarpReduce2(temp_storage2).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
			if_update[v] =1;
			atomicAdd(&gcounter[label-1],1);
			atomicAdd(counter, 1);
			temp_labels[v] = label - 1;
		 }
	 }
}

__global__ void l_tiny_update_syn2_p
( int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int *counter,int* warp_v,int* warp_begin,int warpnumber,double * prob,double* rand_prob)
{
	auto block = cg::this_thread_block();
	auto g = cg::tiled_partition<32>(block);
	const int tid = blockIdx.x*blockDim.x+ threadIdx.x;
	const int warpIdAll = (tid)/ 32;
	int v = -1;
	if(warpIdAll < warpnumber){
		v = warp_v[tid];
	}
	double p = 1-prob[v];
	if(rand_prob[v]<=p){
		v=-1;
	}
	const int begin = warp_begin[warpIdAll];
	// load nodes. group lanes by node id.
	auto mask = __ballot_sync(FULL_MASK, v!=-1);
	int label;
	if(v!= -1){
		label = adj_labels[begin+g.thread_rank()]+1;
	}else{
		label = 0;
	}

	// compute head_flag.
	 auto vmask = __match_any_sync(FULL_MASK,v);
	 auto lmask = __match_any_sync(vmask,label);
	 int count = 0;
	 int head = 0;
	 if(g.thread_rank()==__ffs(vmask)-1){
		 head = 1;
	 }
	 if(g.thread_rank()== __ffs(lmask)-1){
		 count = __popc(lmask);
	 }

	 typedef cub::WarpReduce<unsigned long long> WarpReduce;
	 __shared__ typename WarpReduce::TempStorage temp_storage;
	 unsigned long long valkey = ((unsigned long long)count <<32 )+ label;
	 unsigned long long aggregate = WarpReduce(temp_storage).HeadSegmentedReduce(
			 valkey, head, cub::Max());
	 if(g.thread_rank()==__ffs(vmask)-1 && v!=-1){
		 label =  aggregate;
		 auto lbl = temp_labels[v];
		 if (lbl != label - 1) {
			atomicAdd(counter, 1);
			temp_labels[v] = label - 1;
			prob[v] = 0;
		 }else{
			prob[v] = 0.5+prob[v]/2;
		 }
	 }
}

//template<int VT,int W>
//__global__ void ml_tiny_update_syn3
//(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter, bool* if_update)
//{
//	int mylanekey;
//	int my_max_key;
//	// the vertex about to process
//	#pragma unroll
//	for (int stride = 0;stride < VT;++stride){
//		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)+ offset;
//		if(id< n){
//			float my_max_score = (float)INT_MIN;
//			int v = vertices[id];
//			//printf(" %d \n",id);
//			const int begin = offsets[v];
//			const int end = offsets[v+1];
//			//load value
//			//no neighbor
//			if(begin==end){
//				my_max_key = v+1;
//				temp_labels[v] = my_max_key - 1;
//			}else{
//				// number of neighbors is tiny;
//				int count = 0;
//				float score = 0;
//				for(int j = begin;j<end;j++){
//					mylanekey = adj_labels[j]+1;
//
//					if(score > my_max_score){
//						my_max_score = score;
//						my_max_key = mylanekey;
//					}
//				}
//				auto lbl = temp_labels[v];
//				atomicAdd(&gcounter[my_max_key - 1],1);
//				if (lbl != my_max_key - 1) {
//					if_update[v] = 1;
//					atomicAdd(counter, 1);
//					temp_labels[v] = my_max_key - 1;
//				}
//			}
//		}
//	}
//}
//
//template<int VT>
//__global__ void ml_tiny_update_syn3_cmp
//(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,int* gcounter,int* gcounter_temp, float delta)
//{
//	int mylanekey;
//	int my_max_key;
//
//	// the vertex about to process
//	#pragma unroll
//	for (int stride = 0;stride < VT;++stride){
//		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)+ offset;
//		if(id< n){
//			float my_max_score = (float)INT_MIN;
//			int v = vertices[id];
//			//printf(" %d \n",id);
//			const int begin = offsets[v];
//			const int end = offsets[v+1];
//			//load value
//			//no neighbor
//			if(begin==end){
//				my_max_key = v+1;
//				temp_labels[v] = my_max_key - 1;
//			}else{
//				// number of neighbors is tiny;
//				int count = 0;
//				float score = 0;
//				for(int j = begin;j<end;j++){
//					mylanekey = adj_labels[j]+1;
//					score = count-delta* (float)(gcounter_temp[mylanekey-1]);
//					if(score > my_max_score){
//						my_max_score = score;
//						my_max_key = mylanekey;
//					}
//				}
//
//					temp_labels[v] = my_max_key - 1;
//
//			}
//		}
//	}
//}

//template<int VT>
//__global__ void ml_tiny_update_syn3
//(int offset, int n, int *neighbors,int *offsets, int *adj_labels,int *temp_labels,int* vertices,int *counter,int* gcounter,int* gcounter_temp, float delta)
//{
//	int mylanekey;
//	int my_max_key;
//
//	// the vertex about to process
//	#pragma unroll
//	for (int stride = 0;stride < VT;++stride){
//		auto id = (threadIdx.x + VT*blockIdx.x * blockDim.x + stride*blockDim.x)+ offset;
//		if(id< n){
//			float my_max_score = (float)INT_MIN;
//			int v = vertices[id];
//			//printf(" %d \n",id);
//			const int begin = offsets[v];
//			const int end = offsets[v+1];
//			//load value
//			//no neighbor
//			if(begin==end){
//				my_max_key = v+1;
//				temp_labels[v] = my_max_key - 1;
//			}else{
//				// number of neighbors is tiny;
//				int count = 0;
//				float score = 0;
//				for(int j = begin;j<end;j++){
//					mylanekey = adj_labels[j]+1;
//					score = count-delta* (float)(gcounter_temp[mylanekey-1]);
//					if(score > my_max_score){
//						my_max_score = score;
//						my_max_key = mylanekey;
//					}
//				}
//				auto lbl = temp_labels[v];
//				atomicAdd(&gcounter[my_max_key - 1],1);
//				if (lbl != my_max_key - 1) {
//					atomicAdd(counter, 1);
//					temp_labels[v] = my_max_key - 1;
//				}
//			}
//		}
//	}
//}

// update all values in list
template<typename V, typename E, int VT>
__global__ void ml_gather_all(const V *neighbors, const V *labels, V *labellist, V* gcounter,V* gcounterlist,E m)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];
			auto my_label = labels[my_neighbor];
			labellist[j] =my_label;
			auto count = gcounter[my_neighbor];
			gcounterlist[j] = count;
		}
	}
}
// update all values in list
template<typename V, typename E, int VT>
__global__ void ml_gather_labels(const V *neighbors, const V *labels, V *labellist,E m)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];
			auto my_label = labels[my_neighbor];
			labellist[j] =my_label;

		}
	}
}


// update values in list according to bool array
template<typename V, typename E, int VT>
__global__ void ml_gather_labels(const V *neighbors, const V *labels, V *labellist,E m, bool * if_update)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];
			if(if_update[my_neighbor]){
				auto my_label = labels[my_neighbor];
				labellist[j] =my_label;
			}
		}
	}
}



// update all values in gcounterlist based on labellist.
template<typename V, typename E, int VT>
__global__ void ml_gather_counter(const V *labellist, V* gcounter,V* gcounterlist,E m)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto label = labellist[j];
			auto count = gcounter[label];
			gcounterlist[j] = count;
		}
	}
}

//update according to bool array
template<typename V, typename E, int VT>
__global__ void ml_gather_all2 (const V *neighbors, const V *labels, V *labellist, V* gcounter,V* gcounterlist,V m, bool* if_update,bool* if_update_gcounter)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto my_neighbor = neighbors[j];

			if(if_update[my_neighbor]){
				auto my_label = labels[my_neighbor];
				labellist[j] =my_label;
				auto count = gcounter[my_label];
				gcounterlist[j] = count;
			}else{
				auto my_label = labellist[j];
				if(if_update_gcounter[my_label]){
					auto count = gcounter[my_label];
					gcounterlist[j] = count;
				}
			}
		}
	}
}
template<typename V, typename E, int VT>
__global__ void ml_gather_counter2 ( V* neighbors,V* labellist, V* gcounter,V* gcounterlist,V m, bool* if_update, bool* if_update_gcounter)
{
	#pragma unroll
	for (int i = 0;i< VT;i++){
		int j = threadIdx.x + blockIdx.x * blockDim.x + blockDim.x*gridDim.x*i;
		if(j < m)
		{
			auto neighbor = neighbors[j];
			auto my_label = labellist[j];
			if(if_update[neighbor]){
				auto count = gcounter[my_label];
				gcounterlist[j] = count;
			}else{
				if(if_update_gcounter[my_label]){
					auto count = gcounter[my_label];
					gcounterlist[j] = count;
				}
			}
		}
	}
}



// compute the bool array.
// if gcounter[label] need to update, set if_update_gcounter true.
__global__ void compute_if_gcounter (int* gcounter,int* gcounter_temp,bool* if_update_gcounter,float alpha,int n)
{

		int j = threadIdx.x + blockIdx.x * blockDim.x;
		if(j < n)
		{
			int diff = gcounter[j]-gcounter_temp[j];
			if(diff>alpha || diff<-alpha){
				if_update_gcounter[j] = true;
			}
		}

}

template<int TS>
__global__ void compute_num_blocks // per vertex
(int *offsets, int start, int end, int *num_blocks)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int j = tid +start;
    if(j<end){
    	num_blocks[tid] = (offsets[j + 1] - offsets[j] + TS - 1) / TS;
    }
}


__global__ void assign_blocks
(int *num_blocks, int left, int right, int *block_v, int* block_id)
{
    const int nb = num_blocks[right-left];
    for (auto b: grid_stride_range(nb)) {
        // Binary search
        int d_left = 0;
        int d_right = right-left;
        while (d_left < d_right) {
            int mid = (d_left + d_right) / 2;
            if (num_blocks[mid] <= b) {
                d_left = mid + 1;
            } else {
                d_right = mid;
            }
        }

        block_v[b] =  d_left+ left-1;

        block_id[b] = b - num_blocks[d_left - 1];


    }
}




__global__ void storelabel (int *temp_labels,int* d_label_container, int n,int label_size){
	const int v= get_gid();
	if(v < n ){
//		for (int k = 0; k < label_size; ++k){
//			if(d_label_container[v*label_size] == -1){
				d_label_container[v*label_size] = temp_labels[v];
//				break;
//			}
//		}
	}
}

__global__ void picklabel (int *temp_labels,int* d_label_container, int n,int label_size){
	const int v= get_gid();
	if(v < n ){
		temp_labels[v] = d_label_container[v*label_size];
	}
}
