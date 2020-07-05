// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

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

void writeedge(int u, int v, bool if_txt,std::string filename, std::ofstream& ofs, int* writebuf){
            // write result to file
        if(if_txt == false){
            // write to 32-bit. so number of nodes should no larger than a MAX_INT
            writebuf[0] = u;
            writebuf[1] = v;
            ofs.write((char *) writebuf, sizeof(int) * 2);             

        }else{
            
            unsigned int mark = 0;    
            ofs << u <<' '<< v<<endl;     
        }
}

// build graph: #node, #edge, average degree.
void graphGenerator(std::string filename,unsigned int nodenum ,unsigned int edgenum, bool if_txt,int randomseed){
    std::ofstream ofs;
    if(if_txt == false){
        ofs.open(filename, std::ofstream::out | std::ofstream::binary | std::ofstream::app);
    }else{
         ofs.open(filename, std::ofstream::out | std::ofstream::app);
    }
    //write buff for bin file
    int writebuf[2];
    srand (randomseed);
    int realedge = 0;
    cout << "graph generating ... " <<filename <<endl;
    // p =  m/ n^2 
    double p = (double)edgenum/((double)nodenum *nodenum) *2;
    // cout << "p: "<< p<<endl;

    //classic random graph generator, Edge (u,v)
    int u = 1;
    int v = -1;
    while (u < nodenum)
    {
        /* code */
        double r = rand()/double(RAND_MAX);
        double delta = log(1-r) / log(1-p);
        v = v +1 + (int)delta;
        // cout << "delta: "<<delta<<endl;
        if (v >= u && u < nodenum){
            int time = v/u;
            v = v%u;
            u += time;
        }
        if(u < nodenum){
            // cout <<"u : "<< u <<" v: "<< v<<endl;
            realedge++;
            writeedge (u,v,if_txt,filename,ofs,writebuf);
            if((realedge & 0xffffff) == 0){
                cout << realedge << " / " << edgenum <<endl;
            }
        }
    }
    
    
    // cout << "realedge: " << realedge <<endl;
}



int main(int argc, char *argv[])
{

    if (argc != 6) {
        printf("input argc %d needed %d\n", argc, 6 );

        printf("usage: ./graphGenerator filename #node #edge saveAsTxt randomSeed\n");
        return 0;
    }
    string outputfile(argv[1]);
    unsigned int nodenum = atol(argv[2]);
    unsigned int edgenum = atof(argv[3]);
    bool if_txt = (bool)atoi(argv[4]);
    int randomSeed = atoi(argv[5]);
    graphGenerator(outputfile, nodenum ,edgenum, if_txt, randomSeed);

    return 0;
}
