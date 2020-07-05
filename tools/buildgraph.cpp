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

using namespace std;
// build graph: #node, #edge, average degree.
void graphGenerator(std::string filename,unsigned long long nodenum ,double average_degree, bool if_txt){

    cout << "graph generating ..." <<endl;

    vector<set <unsigned long long>*> edgelists;
    for (auto k = 0;k< nodenum;++k){
         set <unsigned long long >* edge  = new set < unsigned long long >;
        edgelists.push_back(edge);
    }
    double p = average_degree/nodenum;
    unsigned long long edge_num = 0;
    //generate a list of degrees, links to n-1 nodes
    for(unsigned long long  i = 0;i <nodenum-1; ++i){
         unsigned long long u = i;
         unsigned long long v = i+1+rand()%(nodenum - i);
         set <unsigned long long >* edge  = edgelists.at(u);
         // each node at least has one neighbor
         edge->insert(v);
         edge_num++;
        // insert (v,u)
        edgelists.at(v)->insert(u);
         for(unsigned long long  j = i+1;j <nodenum ; ++j){        
            if(j == v){
                continue;
            }
            //because p = average_degree/nodenum
            if(rand()%(nodenum-1) < average_degree-1){
                edge->insert(j);
                edgelists.at(j)->insert(u);
                edge_num++;
            }
         }
    }
    // reorder nodes by its degree.
    cout << "reorder nodes by its degree...";
    std::vector<NodeAndDegree> nd(nodenum);
    for(auto i = 0; i< nodenum; ++ i){
        nd[i].node = i;
        nd[i].degree = edgelists[i]->size();
    }
    sort(nd.begin(), nd.end(),comparator);
    // find the mapping relationship
    std::vector<int> vertex_list(nodenum);
    map<int, int> vertex_mapping;
    for(auto i = 0; i< nodenum; ++ i){
        vertex_list[i] = nd[i].node;
        vertex_mapping[vertex_list[i]] = i;
    }
    //clear memory
    nd.clear();
    nd.shrink_to_fit();

    // write result to file
    if(if_txt == false){
        std::ofstream ofs(filename, std::ofstream::out | std::ofstream::binary);
        int writebuf[2];
        for (auto i = 0; i < nodenum; ++i) {
            if ((i & 0xffffff) == 0) {
                printf("%d / %d\n", i + 1, nodenum);
            }
            auto u = vertex_list[i];
            vector <unsigned long long>edges_of_u (edgelists.at(u)->begin(),edgelists.at(u)->end());
            auto x = vertex_mapping[u];
            writebuf[0] = x;
            for (int j = 0; j < edges_of_u.size(); ++j) {
                auto y = vertex_mapping[edges_of_u[j]];
                if (x >= y) continue;
                writebuf[1] = y;
                ofs.write((char *) writebuf, sizeof(int) * 2);
                
            } 
        }
    }else{

        std::ofstream ofs(filename, std::ofstream::out);

        for (auto i = 0; i < nodenum; ++i) {
            if ((i & 0xffffff) == 0) {
                printf("%d / %d\n", i + 1, nodenum);
            }
            if(i!=0){
                ofs<<std::endl;
            }
            auto u = vertex_list[i];
            vector <unsigned long long>edges_of_u (edgelists.at(u)->begin(),edgelists.at(u)->end());
            auto x = vertex_mapping[u];
           
            for (int j = 0; j < edges_of_u.size(); ++j) {
                auto y = vertex_mapping[edges_of_u[j]];
                if (x >= y) continue;
                
                ofs << x <<' '<< y;
                
            } 

        }

    }
}



int main(int argc, char *argv[])
{

    if (argc != 2) {
        printf("usage: ./graphGenerator  \n");
        return 0;
    }

    const std::string filename("example.txt");

    graphGenerator(filename, 40 ,3, true);


    return 0;
}
