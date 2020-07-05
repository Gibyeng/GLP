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

// DFS 
void DFSutil(int v , map<int,bool> & visited,const std::vector<std::vector<int>> & adj,std::vector <int> &components){
    // mark v as visited 
    visited [v] = true;
    components.push_back(v); 
    for(int i = 0;i< adj[v].size();++i){
        if(!visited[adj[v].at(i)]){
            DFSutil(adj[v].at(i),visited,adj,components);
        }
    }
}

void DFSutil_nonrecursive(int v , map<int,bool> & visited,const std::vector<std::vector<int>> & adj,std::vector <int> &components){
    // mark v as visited 
    visited [v] = true;
    stack<int> my_stack;
    stack<unsigned long> position;
    my_stack.push(v);
    unsigned long i = 0;
    position.push(i);
    components.push_back(v); 
    while (!my_stack.empty()){
        int node = my_stack.top();
        unsigned long start  = position.top();
        for(i = start;i< adj[node].size();++i){
            int cur = adj[node].at(i);
            if(!visited[cur]){
                visited [cur] = true;
                my_stack.push(cur);
                position.pop();
                position.push(i+1);
                 position.push(0);
                components.push_back(cur); 
                break;
            }
        }
        if(i==adj[node].size()){
            my_stack.pop();
            position.pop();
            printf("pop %d \n",node);
        }
    }
}



// Fills Stack with vertices according to DFS
void fillOrder(int v,  map<int,bool> &visited, stack<int> &Stack,const std::vector<std::vector<int>> & adj) 
{ 
    // Mark the current node as visited 
    visited[v] = true; 
  
    // Recur for all the vertices adjacent to this vertex 
    for(int i = 0;i<adj[v].size();++i){
        if(!visited[adj[v].at(i)]) {
            fillOrder(adj[v].at(i), visited, Stack,adj);
        }

    }
    
    // All vertices reachable from v are processed by now, push v  
    Stack.push(v);

}

void fillOrder_nonrecursive(int v,  map<int,bool> &visited, stack<int> &Stack,const std::vector<std::vector<int>> & adj) 
{ 
    // Mark the current node as visited 
    visited[v] = true; 
    // store stack by hand
    stack<int> mystack;
    mystack.push(v);
    Stack.push(v);
    //find one of the vertices adjacent to this vertex 
    while (!mystack.empty()){
        int node = mystack.top();
        unsigned long i = 0;
        for(i = 0;i<adj[node].size();++i){
            int cur = adj[node].at(i);
            if(!visited[cur]) {
                visited[cur] = true; 
                mystack.push(cur);
                Stack.push(cur);
                break;
            }
        }
        if(i==adj[node].size()){
            mystack.pop();
        }     
    }  
}


void find_component(const string &filename){
    std::vector<std::vector<int>> adjacent_list;
    std::vector<std::vector<int>> adjacent_list_reverse;
    std::set<int> vertex_set;
    std::string s;
    long progress = 0;
    std::ifstream ifs(filename);
    while (getline(ifs, s)) {
        if (s[0] == '#' || s[0] == '%') {
            continue;
        }
        int u, v;
        // swap u and v
        sscanf(s.data(), "%d\t%d\n", &v, &u);
        if (v == u) {
            continue;
        }
        // overflow
        if (v<0 || u<0) {
            continue;
        }
        //node set
        vertex_set.insert(u);
        vertex_set.insert(v);
        int max_node = max(u,v);
        if (max_node >= adjacent_list.size()) {
            adjacent_list.resize(max_node + 1);
            
        }
        if(max_node >= adjacent_list_reverse.size()) {
            adjacent_list_reverse.resize(max_node + 1);    
        }
        //direct edge.
        adjacent_list[u].push_back(v);
        //reverse edge
        adjacent_list_reverse[v].push_back(u);

        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    std::vector<int> vertex_list(vertex_set.begin(), vertex_set.end());
    std::cout << "read graph data completed..."<< std::endl;
    stack<int> Stack;
    // if v is visited.
    map<int,bool> visited;
    
    for(int i= 0;i<vertex_list.size();++i){
        visited.insert(pair<int,bool>( vertex_list[i], false) );
    }
    for (int i= 0;i< vertex_list.size();i++){
        if(visited[vertex_list[i]] ==false){
            fillOrder_nonrecursive(i,visited,Stack,adjacent_list);
        }
    }
    cout << "2nd DFS begin ..."<< endl;
    //prepare for the second DFS.
    for(int i= 0;i<vertex_list.size();++i){
        visited[vertex_list[i]] = false;
    }
    // Now process all vertices in order defined by Stack 
    vector<int> components;
    vector<int> max_component;

    int max_c = 0;
    progress = 0;
    while (Stack.empty()==false){
        components.clear();
        int v =Stack.top();
        Stack.pop();
        if(visited[v]== false ){
            DFSutil_nonrecursive(v,visited,adjacent_list_reverse, components);
            int cs = components.size();
            if(max_c < cs){
                max_c = cs;
                max_component = components;
            }
        }
        if (((++progress) & 0xffffff) == 0) {
            printf("%ld\n", progress);
        }
    }
    cout << "remapping begin ..."<< endl;
    // remap node id and sort by degree
    // find the degree of each vertex.
    std::vector<NodeAndDegree> nd(max_component.size());
    for(int i = 0;i < max_component.size();++i){
        nd[i].node = max_component[i];
        nd[i].degree =  adjacent_list[max_component[i]].size();
    }
    //sort by degree
    std::sort(nd.begin(), nd.end(),comparator);
    // copy result back
   for(int i = 0;i < max_component.size();++i){
       max_component[i] = nd[i].node;  
    }
    //clear memory
    nd.clear();
    nd.shrink_to_fit();
    // reordering 
    std::map<int, int> vertex_mapping;
    for (int i = 0; i < max_component.size(); ++i) {
        vertex_mapping[max_component[i]] = i;
    }
    std::string output_name =std::string("sorted_direct_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");

    std::ofstream ofs(output_name, std::ofstream::out | std::ofstream::binary);
    
    int writebuf[2];
    cout << "writing process..."<<endl;
    for (int i = 0; i < max_component.size(); ++i) {
        if ((i & 0xffffff) == 0) {
            printf("%d / %d\n", i + 1, max_component.size());
        }
        auto u = max_component[i];
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
        printf("usage: ./findComponents filename\n");
        return 0;
    }

    const std::string filename(argv[1]);

    find_component(filename);

    return 0;
}
