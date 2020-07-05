// -*- coding: utf-8 -*-
#pragma once
#include "../common/range.cuh"
#include <memory>
# include <cmath>
# include <limits>
#define LONG_PRIME 4294967311


template<int C>
class SharedFA {
public:
    uint32_t keys[C];
    uint32_t vals[C];


    __device__ uint32_t increment(uint32_t key, uint32_t x=1) {
        auto c = _increment(keys, vals, 0, C, key, x);
        return c;
    }

    __device__ void clear() {
        for (auto i: block_stride_range(C)) {
            keys[i] = 0;
            vals[i] = 0;
        }
    }
	__device__ uint32_t hash(uint32_t k) {
		   k ^= k >> 16;
		   k *= 0x85ebca6b;
		   k ^= k >> 13;
		   k *= 0xc2b2ae35;
		   k ^= k >> 16;
		   return k;
	   }
//  first implement
//	__device__ uint32_t _increment
//	    (uint32_t *keys, uint32_t *vals, uint32_t begin, uint32_t end, uint32_t key,uint nt = 64, uint32_t x=1)
//	    {
//	        auto s = end - begin;
//	        auto v = 0;
//	        // if hash function is not full, add a value
//	        if(nitems<C-nt){
//	        	for (uint32_t i = hash(key) % s;; ++i) {
//	            uint32_t k = keys[begin + i];
//					if (k != key) {
//						if (k != 0) {
//							continue;
//						}
//
//						uint32_t prev = atomicCAS(&keys[begin + i], 0, key);
//						if ((prev != 0) && (prev != key)) {
//							continue;
//						}
//					}
//					v = atomicAdd(&vals[begin + i], 1);
//					atomicCAS(&keys[begin + i], 1, key);
//					break;
//				}
//
//	        }else{	// hash function is full,all count -1
//
//	        	for(uint32_t i =begin;i<end;++i){
//	        		auto prev = atomicCAS(&keys[begin + i], keys[begin + i]>0, key-1);
//	        		if(prev == 1){
//	        			atomicCAS(&keys[i], 1, 0);
//	        			atomicAdd(&nitems, -1);
//	        		}
//	        	}
//	        }
//
//	        return v + x;
//
//
//	    }

	__device__ uint32_t _increment (uint32_t *keys, uint32_t *vals, uint32_t begin, uint32_t end, uint32_t key, uint32_t x=1) {
	    	// find the index for insert
	    	 auto s = end - begin;
	    	 auto v = 0;
	    	 bool isfind = 0;
	    	 for (uint32_t i = hash(key), step = 0;step < s; ++i,++step) {
	    		 i = i% s;
				 uint32_t prev = atomicCAS(&keys[begin + i], 0, key);
				 if ((prev != 0) && (prev != key)) {
					continue;
				 }
				 isfind = 1;
				 v = atomicAdd(&vals[begin + i], 1);
				 break;
	    	 }
	    	 //check whether buffer is full, if full all min 1 and reset some value
	    	 if(isfind == 0){
	    		 //printf("full v is %d \n",key);
	    		 for(uint32_t j =begin;j<end;++j){
					auto prev = atomicDec(&vals[j], 0);
					if(prev == 1){
						//printf("j: %d,set to 0,vals : %d,keys: %d \n",j,vals[j],keys[j]);
						atomicExch(&keys[j], 0);
					}
					//printf("j: %d,after vals: %d keys: %d,prev: %d\n",vals[j],keys[j],prev);
	    		 }
	    	 }
	    	 return v + x;
	}
};
