// -*- coding: utf-8 -*-

#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <memory>
#include "myutil.h"



template<typename V, typename E>
class Graph {
    // Virtual class representing a graph
    // * V is the type of vertices
    // * E is the type of edges

public:
    struct iterator {
        iterator(V *first, V *last): _begin(first), _end(last) { }
        iterator(V *first, std::ptrdiff_t size): iterator(first, first + size) {}

        constexpr V *begin() const noexcept { return _begin; }
        constexpr V *end() const noexcept { return _end; }

        V *_begin;
        V *_end;
    };

    Graph() = default;
    Graph(const std::string &filename);
    virtual ~Graph() = default;

    virtual Range<V> iterate_vertices() const = 0;
    virtual iterator iterate_neighbors(V v) = 0;

    V n;  // Number of vertices
    E m;  // Number of undirected edges
    int startvalue; // startswith 1 (1,n) or 0 (0,n-1);
    bool direct; // whether is a direct graph.
};


template<typename V>
using AdjacencyList = std::vector<std::vector<V>>;


// Returns a graph as an adjacency list
template<typename V>
AdjacencyList<V> load_edgelist_txt(const std::string &filename,int startvalue);
template<typename V>
AdjacencyList<V> load_edgelist_bin(const std::string &filename,int startvalue,bool direct, int & max_u_v);
template<typename V>
AdjacencyList<V> load_metisgraph_txt(const std::string &filename,int startvalue);


template<typename V, typename E>
class CSRGraph: public Graph<V, E> {
    // Class representing a graph by the CSR format

public:
    using iterator = typename Graph<V, E>::iterator;

    CSRGraph() = default;
    CSRGraph(const std::string &filename,int startvalue,bool direct) {
    	this->startvalue = startvalue;
    	this->direct = direct;
        auto pos = filename.find_last_of(".");
        auto ext = filename.substr(pos);
        AdjacencyList<V> adj_list;
        if (ext == ".txt" || ext == ".dat"||ext == ".bin") {
        	 std::cout<<"read start"<<std::endl;
        	 V max_u_v = 0;
        	if(ext != ".bin"){

        		adj_list = load_edgelist_txt<V>(filename, this->startvalue);
        	}else{

        		adj_list = load_edgelist_bin<V>(filename, this->startvalue, this->direct, max_u_v);
        	}
        	 std::cout<<"read end"<<std::endl;
        	if(!direct){
        		 this->n = adj_list.size();
        	 }else{
        		 this->n = max_u_v;
        	 }

            auto size_accumulator = [](const E &a, decltype(adj_list[0]) &b) {
                return a + b.size();
            };

            this->m = std::accumulate(adj_list.begin(), adj_list.end(), E(), size_accumulator);

            neighbors.resize(this->m);
            offsets.resize(this->n + 1);
            std::cout << " Adj to CSR begin" <<std::endl;
//            std:: cout << "n "<<this->n <<std::endl;
            // Adjacency list to CSR format
            E cur = 0;
			if(!direct){
				for (auto u: range(this->n)) {
					bool f = false;
					for (auto v: adj_list[u]) {
						if (f && u < v) {
							// // Insert a loop
							// neighbors[cur++] = u;
							f = false;
						}
						neighbors[cur++] = v;
					}
					if (f) {
						neighbors[cur++] = u;
					}
					offsets[u + 1] = cur;
					adj_list[u].clear();
				}
			}else{
				size_t adj_list_size = adj_list.size();
				for (auto u: range(adj_list_size)) {
					for (auto v: adj_list[u]) {
						neighbors[cur++] = v;
					}
					offsets[u + 1] = cur;
					adj_list[u].clear();
				}
			}

//            std::cout <<"offset size: "<<offsets.size();
//            for(int j :range(offsets.size())){
//            	std::cout <<" "<<offsets[j];
//            }
        }
//        else {
//            load_symmetric_edgelist_bin(filename,this->startvalue);
//            this->n = offsets.size() - 1;
//            this->m = neighbors.size();
//        }

//        cudaHostRegister((void *) &neighbors[0], sizeof(V) * this->m,
//                         cudaHostRegisterMapped | cudaHostRegisterPortable);
//        cudaHostRegister((void *) &offsets[0], sizeof(E) * (this->n + 1),
//                         cudaHostRegisterMapped | cudaHostRegisterPortable);
    }

    ~CSRGraph() {
//        cudaHostUnregister((void *) &neighbors[0]);
//        cudaHostUnregister((void *) &offsets[0]);
    }

    void load_symmetric_edgelist_bin(const std::string &filename,int startvalue) {
        std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
        ifs.seekg(0, std::ios::beg);
        const int bufsize = 2;
        V buf[bufsize];

        int cur = -1;
        bool f = false;
        std::cout<<"read start"<<std::endl;
        while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
            V u = buf[0]-startvalue;
            V v = buf[1]-startvalue;
            //printf("U: %d,V: %d \n",u,v);
            if (u != cur) {

                offsets.push_back(neighbors.size());
                cur = u;
                f = true;
            }
            if (f && u < v) {
                // // Insert a loop
                // neighbors.push_back(u);
                f = false;
            }

            neighbors.push_back(v);
        }

        offsets.push_back(neighbors.size());
        ifs.close();
        std::cout<<"read completed"<<std::endl;
    }

    Range<V> iterate_vertices() const {
        return range(0, this->n);
    }

    iterator iterate_neighbors(V v) {
        return iterator(&neighbors[offsets[v]], &neighbors[offsets[v + 1]]);
    }

    std::vector<V> neighbors;
    std::vector<E> offsets;

};



template<typename V>
AdjacencyList<V> load_edgelist_txt(const std::string &filename,int startvalue)
{
	  std::ifstream ifs(filename);
	    AdjacencyList<V> adj_list;
	    // Read the file line by line
	    std::string s;
	    while (getline(ifs, s)) {
	        V u, v;
	        std::stringstream ss(s);
	        ss >> u >> v;
	        u = u - startvalue;
	        v = v - startvalue;
	        if (adj_list.size() <= v) {
				adj_list.resize(v + 1);
			}
        	adj_list[u].push_back(v);
			// Symmetrize
			// * Assume the edge (u, v), where u > v, is not in the file
			adj_list[v].push_back(u);

    }
    return adj_list;
}
template<typename V>
AdjacencyList<V> load_edgelist_bin(const std::string &filename,int startvalue,bool direct, int & max_u_v)
{
	std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
	V buf[bufsize];
	 AdjacencyList<V> adj_list;
    // Read the file line by line
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
    	V u = buf[0]-startvalue;
		V v = buf[1]-startvalue;
		if(max_u_v < u|| max_u_v < v){
			max_u_v = max(u,v);
		}
        if(u!=v ){
        	if (adj_list.size() <= u) {
				 adj_list.resize(u + 1);
			 }
			adj_list[u].push_back(v);

			if(!direct){
				if (adj_list.size() <= v) {
					 adj_list.resize(v + 1);
				}
				adj_list[v].push_back(u);
			}
        }

    }

    return adj_list;
}


