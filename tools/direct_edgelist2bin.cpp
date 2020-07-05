// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <list> 
#include <stack> 
using namespace std;
// used for sorting
struct NodeAndDegree{
    int node;
    int degree;
};
struct Comparator {
    bool operator () (NodeAndDegree const &i,NodeAndDegree const &j) {
        
        return i.degree > j.degree;
    }
}comparator;


void buildfromtxt(const string &filename){
    std::vector<std::vector<int>> adjacent_list;
    std::set<int> vertex_set;
    std::string s;
    long progress = 0;
    std::ifstream ifs(filename);
    while (getline(ifs, s)) {
        if (s[0] == '#'|| s[0] == '%') {
            continue;
        }
        int u, v;

        //sscanf(s.data(), "%d\t%d\n", &v, &u);
        sscanf(s.data(), "%d\t%d\n", &u, &v);
        std::cout<< "u: "<<u << " v: "<< v<<std::endl;
        if (v == u) {
            continue;
        }
        //node set
        if( u < 5013066 || v < 5013066){
            vertex_set.insert(u);
            vertex_set.insert(v);
            int max_node = max(u,v);
            if (max_node >= adjacent_list.size()) {
                adjacent_list.resize(max_node + 1);
            
        }
        
        //direct edge.
        adjacent_list[u].push_back(v);
        }
       

        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
    std::cout << "read graph data completed..."<< std::endl;


    cout << "remapping begin ..."<< endl;
    // remap node id and sort by degree
    // find the degree of each vertex.
    std::vector<NodeAndDegree> nd(vertex_list.size());
    for(int i = 0;i < vertex_list.size();++i){
        nd[i].node = vertex_list[i];
        nd[i].degree =  adjacent_list[vertex_list[i]].size();
    }
    //sort by degree
    std::sort(nd.begin(), nd.end(),comparator);
    // copy result back
   for(int i = 0;i < vertex_list.size();++i){
       vertex_list[i] = nd[i].node;  
    }
    //clear memory
    nd.clear();
    nd.shrink_to_fit();
    // reordering 
    cout << "reordering begin ..."<< endl;
    std::map<int, int> vertex_mapping;
    for (int i = 0; i < vertex_list.size(); ++i) {
        vertex_mapping[vertex_list[i]] = i;
    }
    std::string output_name =std::string("sorted_direct_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    int writebuf[2];
    cout << "writing process..."<<endl;
    for (int i = 0; i < vertex_list.size(); ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, vertex_list.size());
        }
        auto u = vertex_list[i];
        std::sort(adjacent_list[u].begin(), adjacent_list[u].end());
        auto x = vertex_mapping[u];
        writebuf[0] = x;
        for (int j = 0; j < adjacent_list[u].size(); ++j) {
            if(vertex_mapping.find(adjacent_list[u][j])==vertex_mapping.end()){
                
                vertex_mapping[adjacent_list[u][j]] = vertex_mapping.size();
            }
            auto y = vertex_mapping[adjacent_list[u][j]];
            writebuf[1] = y;
            ofs.write((char *) writebuf, sizeof(int) * 2); 
        } 
    }

}

void buildfrombin(const string &filename){
    std::vector<std::vector<int>> adjacent_list;
    std::set<int> vertex_set;
    std::string s;
    long progress = 0;
    const int bufsize = 2;
    int buf[bufsize];
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        int u, v;
        u = buf[0];
		v = buf[1];
        //sscanf(s.data(), "%d\t%d\n", &v, &u);
        sscanf(s.data(), "%d\t%d\n", &u, &v);

        if (v == u) {
            continue;
        }
        //node set
        if( u < 5013066 || v < 5013066){
            vertex_set.insert(u);
            vertex_set.insert(v);
            int max_node = max(u,v);
            if (max_node >= adjacent_list.size()) {
                adjacent_list.resize(max_node + 1);
            
        }
        
        //direct edge.
        adjacent_list[u].push_back(v);
        }
       

        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
    std::cout << "read graph data completed..."<< std::endl;


    cout << "remapping begin ..."<< endl;
    // remap node id and sort by degree
    // find the degree of each vertex.
    std::vector<NodeAndDegree> nd(vertex_list.size());
    for(int i = 0;i < vertex_list.size();++i){
        nd[i].node = vertex_list[i];
        nd[i].degree =  adjacent_list[vertex_list[i]].size();
    }
    //sort by degree
    std::sort(nd.begin(), nd.end(),comparator);
    // copy result back
   for(int i = 0;i < vertex_list.size();++i){
       vertex_list[i] = nd[i].node;  
    }
    //clear memory
    nd.clear();
    nd.shrink_to_fit();
    // reordering 
    cout << "reordering begin ..."<< endl;
    std::map<int, int> vertex_mapping;
    for (int i = 0; i < vertex_list.size(); ++i) {
        vertex_mapping[vertex_list[i]] = i;
    }
    std::string output_name =std::string("sorted_direct_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    int writebuf[2];
    cout << "writing process..."<<endl;
    for (int i = 0; i < vertex_list.size(); ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, vertex_list.size());
        }
        auto u = vertex_list[i];
        std::sort(adjacent_list[u].begin(), adjacent_list[u].end());
        auto x = vertex_mapping[u];
        writebuf[0] = x;
        for (int j = 0; j < adjacent_list[u].size(); ++j) {
            if(vertex_mapping.find(adjacent_list[u][j])==vertex_mapping.end()){
                
                vertex_mapping[adjacent_list[u][j]] = vertex_mapping.size();
            }
            auto y = vertex_mapping[adjacent_list[u][j]];
            writebuf[1] = y;
            ofs.write((char *) writebuf, sizeof(int) * 2); 
        } 
    }

}


int main(int argc, char *argv[])
{

    if (argc != 2) {
        printf("usage: ./direct_edgelist2bin filename\n");
        return 0;
    }
   
    const std::string filename(argv[1]);

    buildfrombin(filename);

    return 0;
}
