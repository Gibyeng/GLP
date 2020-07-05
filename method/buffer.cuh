// -*- coding: utf-8 -*-
#pragma once
#include "../common/range.cuh"
#include <memory>
# include <cmath>
# include <limits>
#define LONG_PRIME 4294967311

// C = 32*warp number in a block ,C: length of buffer in one block,
template<int C>
class Buffer {
public:
	 uint32_t keys[C];
	 uint32_t vals[C];

    __device__ uint32_t hash(uint32_t k) {
		   k ^= k >> 16;
		   k *= 0x85ebca6b;
		   k ^= k >> 13;
		   k *= 0xc2b2ae35;
		   k ^= k >> 16;
		   return k;
	   }
    __device__ void clear() {
         for (auto i: block_stride_range(C)) {
             keys[i] = 0;
             vals[i] = 0;
         }
     }
    __device__ uint32_t increment(uint32_t key,uint32_t begin, uint32_t end, uint32_t x=1) {
           auto c = _increment(this->keys, this->vals, begin, end, key, x);

           return c;
       }

    __device__ uint32_t _increment (uint32_t *keys, uint32_t *vals, uint32_t begin, uint32_t end, uint32_t key, uint32_t x=1) {
    	// find the index for insert
    	 auto s = end - begin;
    	 auto v = 0;
    	 bool isfind = 0;
    	 for (uint32_t i = hash(key) , step = 0;step < s; ++i,++step) {
    		 i = i%s;

			 uint32_t k = keys[begin + i];
			 uint32_t prev = atomicCAS(&keys[begin + i], 0, key);
			 if ((prev != 0) && (prev != key)) {
		  		continue;
		  	 }

			 isfind = 1;
			 v = atomicAdd(&vals[begin + i], x);

			 break;
    	 }
    	 //check whether buffer is full, if is full all min x ect.1 and reset some value, not work when one warp one node.
    	 if(isfind == 0){
//    		 printf("all minus 1 cause of %d \n",key-1 );
//    		 printf("after -1 hash table is: \n");
			 for(uint32_t j =begin;j<end;++j){
				 auto prev = atomicExch(&vals[j], vals[j]-x>0?vals[j]-x:0);
				 //auto prev = atomicDec(&vals[j], 0);
				if(prev <=x){
					atomicExch(&keys[j], 0);
				}
//				printf("i: %d,key: %d,val: %d \n",j-begin,keys[j],vals[j]);
			}

		 }
    	 return v + x;

    }
};
