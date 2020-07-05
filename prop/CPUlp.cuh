// -*- coding: utf-8 -*-
#pragma once

#include "../graph.h"
#include "../method/hashtable.cuh"
#include "../method/countmin.cuh"
#include "../kernel.cuh"
#include <omp.h>

template<typename V, typename E>
class CPUlp {
public:
	using GraphT = CSRGraph<V, E>;
	CPUlp(GraphT* G):G(G){};

    double sync_run_1(int niter);
    double layered_run(int niter);
    double slpa_run(int niter);
    void sync_run_2(int niter);
protected:
    // Attributes
    GraphT* G;
    V* labels;

};

//CPU OPENMP LP baseline
template<typename V, typename E>
double CPUlp<V, E>::sync_run_1(int niter){
	//various of CPU classic method
	std::vector<V> vertices(G->n);
	this->labels = new V[G->n];
	auto prelabels = new V[G->n];
	 for (auto i: range(G->n)) {
		vertices[i] = i;
		labels[i] = i;
		prelabels[i] = i;
	}
	 printf("sync, deterministic \n");
	 Timer t1;
	 float t2 = 0.0;
	 for (auto i: range(niter)) {
		 t1.start();
	        V nupdates = 0;
	        V jumpdates = 0;
	        // run cpu
	        for (auto k: range(G->n)) {
	        		prelabels[k] = labels[k];
	        	}
			#pragma omp parallel for num_threads(20)
	        for (int k = 0;k<G->n;k++) {
//	        	std::cout << omp_in_parallel()<<std::endl;
	        	int u = vertices[k];
	            std::map<V, int> label_count;
	            V max_label = prelabels[u];
	            int max_score = 0;

				int begin = G->offsets[u];
				int end = G->offsets[u+1];
				if(end >0){
					for (int index = begin; index <end; ++index) {
						int v = G->neighbors[index];
						//label starts with 0 but graph id starts with 1.
						V label = labels[v];
						++label_count[label];

						int c = label_count[label];
						if(c> max_score){
							max_score = c;
							max_label = label;
						}
					}
				}
	            if (prelabels[u] != max_label) {
	                labels[u] = max_label;
	                ++nupdates;
	            }

	            label_count.clear();

	        }
	        t1.stop();
	        t2 += t1.elapsed_time();
	        printf("Iteration: %d, changed label(1): %d, t1: %f \n, ",i,nupdates,t1.elapsed_time());

			std::cout<<std::endl;
	    }
	 std::cout<<"all time: "<<t2<<std::endl;
	 delete[] labels;
	 delete[] prelabels;
	 return t2;
}

template<typename V, typename E>
void CPUlp<V, E>::sync_run_2(int niter){
	//various of CPU classic method
	std::vector<V> vertices(G->n);
	this->labels = new V[G->n];
	auto prelabels = new V[G->n];
	 for (auto i: range(G->n)) {
		vertices[i] = i;
		labels[i] = i;
		prelabels[i] = i;
	}

	 srand((unsigned)time(NULL));
	 printf("sync, random\n");
	 Timer t1;
	 float t2 = 0.0;
	 for (auto i: range(niter)) {
		 	t1.start();
	        V nupdates = 0;
	        V jumpdates = 0;
	        // run cpu
	        for (auto k: range(G->n)) {

	        		prelabels[k] = labels[k];
	        	}
	        for (auto u: vertices) {

	            std::map<V, int> label_count;
	            V max_label = prelabels[u];
	            int max_score = 0;
	            std::vector<int> candidate;
	            for (auto v: G->iterate_neighbors(u)) {
					//label starts with 0 but graph id starts with 1.
					V label = labels[v];
					++label_count[label];

					int c = label_count[label];
					if(c == max_score){
						//if two label has the same count, randomly pick one as the max_label
						candidate.push_back(label);
						int index = rand()%(candidate.size());
						max_label = candidate[index];
					}
					if(c> max_score){
						// find a new max_label.
						max_score = c;
						max_label = label;
						candidate.clear();
						candidate.push_back(label);
					}

				}

	            if (prelabels[u] != max_label) {
	                labels[u] = max_label;
	                ++nupdates;
	            }


	            label_count.clear();

	        }

	        t1.stop();
			t2 += t1.elapsed_time();
			printf("Iteration: %d, changed label(1): %d, t1: %f \n, ",i,nupdates,t1.elapsed_time());

			std::cout<<std::endl;
	    }
	 std::cout<<"all time: "<<t2<<std::endl;
	 delete[] labels;
	 delete[] prelabels;

}

template<typename V, typename E>
double CPUlp<V, E>::layered_run(int niter){
	int d_iter = 10;
	int delta = 1;
	//various of CPU classic method
	std::vector<V> vertices(G->n);
	this->labels = new V[G->n];
	auto prelabels = new V[G->n];
	auto gcount_r = new V [G->n];
	auto gcount_w = new V [G->n];
	memset(gcount_r, 0, sizeof(int)*G->n);
	memset(gcount_w, 0, sizeof(int)*G->n);
	 for (auto i: range(G->n)) {
		vertices[i] = i;
		labels[i] = i;
		prelabels[i] = i;
	}
	 printf("layered , sync, deterministic \n");
	 Timer t1;
	 float t2 = 0.0;
	 for(int j = 0; j<d_iter; j++){
		 if(j!=0){
			 for (auto i: range(G->n)) {
				 labels[i] = i;
			 }
		 }
		 for (auto i: range(niter)) {
			 t1.start();
				V nupdates = 0;
				V jumpdates = 0;
				delta = 2*delta;
				int delta_star = delta+1;
				// run cpu
				for (auto k: range(G->n)) {
						prelabels[k] = labels[k];
					}
				memcpy(gcount_r,gcount_w, sizeof(int)*G->n);
				memset(gcount_w, 0, sizeof(int)*G->n);

				#pragma omp parallel for num_threads(20)
				for (int k = 0;k<G->n;k++) {
//		        	std::cout << omp_in_parallel()<<std::endl;

					int u = vertices[k];
					std::map<V, int> label_count;
					V max_label = prelabels[u];
					int max_score = INT_MIN;

					int begin = G->offsets[u];
					int end = G->offsets[u+1];


					for (int index = begin; index <end; index++) {
						//label starts with 0 but graph id starts with 1.

						int v = G->neighbors[index];

						V label = labels[v];

						++label_count[label];

						int c = label_count[label];
						int s = c*delta_star;
						if(label < G->n) {
							s = s- gcount_r[label];
						}
						if(s> max_score){
							max_score = s;
							max_label = label;
						}
					}
					if(max_label<G->n){
						gcount_w[max_label]++;
					}
					if (prelabels[u] != max_label) {
						labels[u] = max_label;
						++nupdates;
					}

					label_count.clear();

				}
				t1.stop();
				t2 += t1.elapsed_time();
				printf("Iteration: %d, changed label(1): %d, t1: %f \n, ",i,nupdates,t1.elapsed_time());

				std::cout<<std::endl;
			}
	 }
	 std::cout<<"all time: "<<t2<<std::endl;
	 delete[] labels;
	 delete[] prelabels;
	 return t2;
}

template<typename V, typename E>
double CPUlp<V, E>::slpa_run(int niter){
	//various of CPU classic method
	std::vector<V> vertices(G->n);
	this->labels = new V[G->n];
	auto prelabels = new V[G->n];
	int* label_candidate;
	int label_size = 5;
	label_candidate = new V[G->n*label_size];
	 for (auto i: range(G->n)) {
		vertices[i] = i;
		labels[i] = i;
		prelabels[i] = i;
	}
	#pragma omp parallel for num_threads(20)
	 for (int i = 0;i< G->n*label_size;i++) {
		 label_candidate[i] = -1;
	 }
	 printf("sync, deterministic \n");
	 Timer t1;
	 float t2 = 0.0;
	 for (auto i: range(niter)) {
		 t1.start();
	        V nupdates = 0;
	        V jumpdates = 0;
	        // run cpu
	        //save labels to candidates
		    #pragma omp parallel for num_threads(20)
	        for (int k = 0; k< G->n; k++){
//	        	int pos =0;
//	        	for(int p = 0;p<label_size;++p){
//	        		if(label_candidate[k*label_size + pos + p] == -1){
//	        			label_candidate[k*label_size+pos+p]= labels[k];
//	        			break;
//	        		}
//	        	}
	        	label_candidate[k*label_size]= labels[k];
			}
	        //read labels from candidates
			#pragma omp parallel for num_threads(20)
	        for (int k = 0; k< G->n; k++){
//	        	int pos =rand();
//				for(int p = 0;p<label_size;++p){
//					pos = (pos+p)%label_size;
//					if(label_candidate[k*label_size + pos] == -1){
//						prelabels[k] = label_candidate[k*label_size + p];
//						break;
//					}
//				}
	        	int pos = 0;
	        	prelabels[k] = label_candidate[k*label_size + pos];
			}
			#pragma omp parallel for num_threads(20)
	        for (int k = 0;k<G->n;k++) {
//	        	std::cout << omp_in_parallel()<<std::endl;
	        	int u = vertices[k];
	            std::map<V, int> label_count;
	            V max_label = prelabels[u];
	            int max_score = 0;
	            for (auto v: G->iterate_neighbors(u)) {
	                //label starts with 0 but graph id starts with 1.
	                V label = labels[v];
	                ++label_count[label];

					int c = label_count[label];
					if(c> max_score){
						max_score = c;
						max_label = label;
					}
	            }

	            if (prelabels[u] != max_label) {
	                labels[u] = max_label;
	                ++nupdates;
	            }

	            label_count.clear();

	        }
	        t1.stop();
	        t2 += t1.elapsed_time();
	        printf("Iteration: %d, changed label(1): %d, t1: %f \n, ",i,nupdates,t1.elapsed_time());

			std::cout<<std::endl;
	    }
	 std::cout<<"all time: "<<t2<<std::endl;
	 delete[] labels;
	 delete[] prelabels;
	 return t2;
}
