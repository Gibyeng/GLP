// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"

int cpu_hash1(int k) {
        k ^= k >> 10;
        k *= 0x82ecba63;
        k ^= k >> 12;
        k *= 0xc3b15e32;
        k ^= k >> 9;
        return abs(k);
    }

int  cpu_hash2(int k) {
        k ^= k >> 12;
        k *= 0xe85cb46c;
        k ^= k >> 17;
        k *= 0xc3b2cf55;
        k ^= k >> 12;
        return abs(k);
    }

int  cpu_hash3(int k) {
        k ^= k >> 16;
        k *= 0x85ebca6b;
        k ^= k >> 13;
        k *= 0xc2b2ae35;
        k ^= k >> 16;
        return abs(k);
    }

template<typename V, typename E>
class SMmlp {
public:
	using GraphT = CSRGraph<V, E>;
	SMmlp(GraphT* G):G(G){};
    void run(int niter);
protected:
    // Attributes
    GraphT* G;
    V* labels;

};


// computer error ratio per niter
template<typename V, typename E>
void SMmlp<V, E>::run(int niter){
	//various of CPU classic method
	std::vector<V> vertices(G->n);
	int* g_counter = new V [G->n];
	int* g_counter_temp = new V [G->n];
	int* labels_jump = new V [G->n];
	this->labels = new V[G->n];
	auto prelabels = new V[G->n];
	 for (auto i: range(G->n)) {
		vertices[i] = i;
		labels[i] = i;
		prelabels[i] = i;
		g_counter[i] = 0;
		g_counter_temp[i] = 0;
		labels_jump[i] = i;
	}
	 float delta = 0.5;
	 int big =0;
	 int medium = 0;
	 int small = 0;
	 // compute how many nodes
	 for (auto u: vertices) {
		 if(G->offsets[u+1] -G->offsets[u]<=5*32){
			if(G->offsets[u+1] -G->offsets[u]<=32){
				small++;
			}else{
				medium++;
			}
		}else{
			big++;
		}
	}
	 printf("small: %d, medium: %d, big: %d \n", small, medium, big);
	for (auto i: range(niter)) {

	        V nupdates = 0;
	        V jumpdates = 0;
	        // run cpu
	        for (auto k: range(G->n)) {
	        	    labels_jump[k] = prelabels[k];
	        		prelabels[k] = labels[k];
	            	g_counter_temp[k] = g_counter[k];
	           		g_counter[k] = 0;
	        	}
	        for (auto u: vertices) {
	        	int *cm_1 = new int [256];
	        	int *cm_2 = new int [256];
	        	int *cm_3 = new int [256];
	        	for(int k = 0;k< 256;k++){
	        		cm_1[k] = 0;
	        		cm_2[k] = 0;
	        		cm_3[k] = 0;
	        	}
	            std::map<V, int> label_count;
	            V max_label = prelabels[u];
	            float max_score = (float)INT_MIN;
	            bool filter = false;
	            for (auto v: G->iterate_neighbors(u)) {
	                //label starts with 0 but graph id starts with 1.
	                V label = prelabels[v];
	                float c ;
	                if(G->offsets[u+1] -G->offsets[u]<=5*32){
	                	if(G->offsets[u+1] -G->offsets[u]<=32){
	                		++label_count[label];
	                		c = label_count[label]- delta* (float)(g_counter_temp[label]);

	                	}else{
						 ++label_count[label];
						 c = label_count[label]- delta* (float)(g_counter_temp[label]);

	                	}

	                }else{
	                	int hash1 = cpu_hash1(label+1)%256;
	                	int hash2 = cpu_hash2(label+1)%256;
	                	int hash3 = cpu_hash3(label+1)%256;
	                	int count1 = cm_1[hash1];
	                	cm_1[hash1]++;
	                	int count2 = cm_2[hash2];
	                	cm_2[hash2]++;
	                	int count3 = cm_3[hash3];
	                	cm_3[hash3]++;
	                	int count = min(count1,count2);
	                	count = min(count,count3)+1;
	                	c = (float)count - delta* (float)g_counter_temp[label];
	                }

					if(c> max_score){
						max_score = c;
						max_label = label;
					}
	            }

				if (prelabels[u] != max_label) {
					labels[u] = max_label;
					++nupdates;
				}
				if (labels_jump[u] != max_label) {
					++jumpdates;
				}
				g_counter[max_label] +=1;



	            label_count.clear();
	            delete []cm_1;
	            delete []cm_2;
	            delete []cm_3;

	        }

	        printf("Iteration: %d, changed label(1): %d, changed label(2): %d  \n, ",i,nupdates,jumpdates);
//	        printf(" label list ");
//	        for (int j =0;j< 20;j++){
//
//		        printf(" %d ",labels[j]);
////	        	 printf(" %d ",prelabels[G->neighbors[j]]);
//	        }
//	        std::cout<<std::endl;
//	        printf(" global label list: ");
//	        for (int j =0;j< 20;j++){
//
//				 printf(" %d ",g_counter_temp[j]);
//			}
//
//			std::cout<<std::endl;
	    }
	 delete[] labels;
	 delete[] prelabels;

}
