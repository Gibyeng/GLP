// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>

struct NodeAndDegree{
    int node;
    int degree;
};
struct Comparator {
    bool operator () (NodeAndDegree const &i,NodeAndDegree const &j) {
        
        return i.degree > j.degree;
    }
}comparator;

void normalize(const std::string &filename)
{
    std::vector<std::vector<int>> adjacent_list;
    std::set<int> vertex_set;

    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    int buf[bufsize];
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        if (buf[0] == '#') {
            continue;
        }
        int u, v;
        u = buf[0];
		v = buf[1];
        vertex_set.insert(u);
        vertex_set.insert(v);
        if (v == u) {
            continue;
        }
        if(u>v){
            std::swap(u,v);
        }
        // so u < v
        if (v >= adjacent_list.size()) {
            adjacent_list.resize(v + 1);
        }

        adjacent_list[u].push_back(v);
        adjacent_list[v].push_back(u);
        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }

    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
    std::sort(vertex_list.begin(), vertex_list.end());
    std::map<int, int> vertex_mapping;
    for (int i = 0; i < vertex_list.size(); ++i) {
        vertex_mapping[vertex_list[i]] = i;
    }
    std::string output_name = filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    int writebuf[2];
    for (int i = 0; i < vertex_list.size(); ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, vertex_list.size());
        }
        auto u = vertex_list[i];
        std::sort(adjacent_list[u].begin(), adjacent_list[u].end());
        auto x = vertex_mapping[u];
        writebuf[0] = x;
        for (int j = 0; j < adjacent_list[u].size(); ++j) {
            auto y = vertex_mapping[adjacent_list[u][j]];
            if (x >= y) continue;
            writebuf[1] = y;
            ofs.write((char *) writebuf, sizeof(int) * 2);
            //std::cout << "U: "<< x << " V: "<< y<<std::endl;
        }
    }
}

// sort graph by degree, remove same edges
void normalize_remove_sort(const std::string &filename)
{
    std::vector<std::vector<int>> adjacent_list;
    std::set<int> vertex_set;

    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    int buf[bufsize];
    int max_nodeid = 0;
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        if (buf[0] == '#') {
            continue;
        }
        int u, v;
        u = buf[0];
		v = buf[1];
        vertex_set.insert(u);
        vertex_set.insert(v);
        if (v == u) {
            continue;
        }
        if(u>v){
            std::swap(u,v);
        }
        // so u < v
        if(v > max_nodeid){
            max_nodeid = v;
        }
        if (v >= adjacent_list.size()) {
            adjacent_list.resize(v + 1);
        }

        adjacent_list[u].push_back(v);
        adjacent_list[v].push_back(u);
        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    std::cout << "node num: " << max_nodeid +1<< "edge num: " << progress;
    // find the degree of each vertex.
    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
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
    std::map<int, int> vertex_mapping;
    for (int i = 0; i < vertex_list.size(); ++i) {
        vertex_mapping[vertex_list[i]] = i;
    }
    std::string output_name =std::string("sorted_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    int writebuf[2];
    for (int i = 0; i < vertex_list.size(); ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, vertex_list.size());
        }
        auto u = vertex_list[i];
        std::sort(adjacent_list[u].begin(), adjacent_list[u].end());
        auto x = vertex_mapping[u];
        writebuf[0] = x;
        for (int j = 0; j < adjacent_list[u].size(); ++j) {
            auto y = vertex_mapping[adjacent_list[u][j]];
            if (x >= y) continue;
            writebuf[1] = y;
            ofs.write((char *) writebuf, sizeof(int) * 2);
            
        } 
    }
}


// sort graph by degree
void normalize_sort(const std::string &filename)
{
    std::vector<std::vector<int>> adjacent_list;
    std::set<int> vertex_set;

    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    int buf[bufsize];
    int max_nodeid = 0;
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        if (buf[0] == '#') {
            continue;
        }
        int u, v;
        u = buf[0];
		v = buf[1];
        vertex_set.insert(u);
        vertex_set.insert(v);
        if (v == u) {
            continue;
        }
        if(u>v){
            std::swap(u,v);
        }
        // so u < v
         if(v > max_nodeid){
            max_nodeid = v;
        }
        if (v >= adjacent_list.size()) {
            adjacent_list.resize(v + 1);
        }
        

        adjacent_list[u].push_back(v);
        adjacent_list[v].push_back(u);
        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    std::cout << "node num: " << max_nodeid +1<< "edge num: " << progress;
    // find the degree of each vertex.
    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
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
    std::map<int, int> vertex_mapping;
    for (int i = 0; i < vertex_list.size(); ++i) {
        vertex_mapping[vertex_list[i]] = i;
    }
    std::string output_name =std::string("sorted_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    int writebuf[2];
    for (int i = 0; i < vertex_list.size(); ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, vertex_list.size());
        }
        auto u = vertex_list[i];
        std::sort(adjacent_list[u].begin(), adjacent_list[u].end());
        auto x = vertex_mapping[u];
        writebuf[0] = x;
        for (int j = 0; j < adjacent_list[u].size(); ++j) {
            auto y = vertex_mapping[adjacent_list[u][j]];
            if (x >= y) continue;
            writebuf[1] = y;

            
            ofs.write((char *) writebuf, sizeof(int) * 2);
        } 
    }
}

// sort graph by degree save mode
void normalize_sort_save(const std::string &filename)
{
    std::vector<std::vector<int>> adjacent_list;
    std::set<int> vertex_set;
    std::cout << "begin to sort graph..."<<std::endl;
    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    int buf[bufsize];
    int max_nodeid = 0;
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        if (buf[0] == '#') {
            continue;
        }
        int u, v;
        u =abs(buf[0]);
		v =abs(buf[1]);
        // std::cout << "u: " << u << "v: "<<v <<std::endl; 
        vertex_set.insert(u);
        vertex_set.insert(v);
        if (v == u) {
            continue;
        }
        if(u>v){
            std::swap(u,v);
        }
        // so u < v
        if(v> max_nodeid){
            max_nodeid = v;
        }
        if (v >= adjacent_list.size()) {
            adjacent_list.resize(v + 1);
        }

        adjacent_list[u].push_back(v);
        adjacent_list[v].push_back(u);
        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    std::cout << "node num: " << max_nodeid +1<< "edge num: " << progress<<std::endl;
    int vertex_num = vertex_set.size();
    std::cout << "after remove 0-degree vertices, new vertex_num : "<< vertex_num<<std::endl;
    // find the degree of each vertex.
    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
    std::vector<NodeAndDegree> nd(vertex_num);
    for(size_t i = 0;i <vertex_num;++i){
        nd[i].node = vertex_list[i];
        nd[i].degree =  adjacent_list[vertex_list[i]].size();
    }
    //sort by degree
    std::sort(nd.begin(), nd.end(),comparator);
    // copy result back
   for(int i = 0;i < vertex_num;++i){
       vertex_list[i] = nd[i].node;  
    }
    //clear memory
    nd.clear();
    nd.shrink_to_fit();
    // reordering 
    std::map<int , int> vertex_mapping;
    for (int i = 0; i < vertex_num; ++i) {
        vertex_mapping[vertex_list[i]] = i;
    }
    std::string output_name =std::string("sorted_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    int writebuf[2];
    for (size_t i = 0; i < vertex_num; ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, vertex_num);
        }
        auto u = vertex_list[i];
        // std::sort(adjacent_list[u].begin(), adjacent_list[u].end());
        auto x = vertex_mapping[u];
        writebuf[0] = x;
        for (int j = 0; j < adjacent_list[u].size(); ++j) {
            auto y = vertex_mapping[adjacent_list[u][j]];
            if (x >= y) continue;
            writebuf[1] = y;

            // write 32-bit
            ofs.write((char *) writebuf, sizeof(int) * 2);
        } 
    }
}

int main(int argc, char *argv[])
{

    if (argc != 2) {
        printf("usage: ./normalize filename\n");
        return 0;
    }

    const std::string filename(argv[1]);

    // normalize(filename);
    normalize_sort_save(filename);


    return 0;
}
