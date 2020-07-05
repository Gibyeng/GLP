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


using namespace std;

void writeedge(int u, int v, std::ofstream& ofs){
            // write result to file  
            ofs << u <<' '<< v;    
}

// merge graph: read graph, write graph, copy or merge, increase nodeid,
void mergegraph(std::string readfilename,std::string writefilename, int if_copy, int increase, int if_delete){
    cout << "start to merge " <<readfilename << " into " << writefilename<< "..."<<endl;
    std::ifstream ifs;
    std::ofstream ofs;
    ifs.open(readfilename, std::ofstream::in);
    ofs.open(writefilename, std::ofstream::out | std::ofstream::app);
    if(if_copy == 1){
        ofs << ifs.rdbuf();
        ifs.close();
        ofs.close();
    }else{
        std::string s;
        while (getline(ifs, s)) {
	        int u, v;
	        std::stringstream ss(s);
	        ss >> u >> v;
	        u = u + increase;
	        v = v + increase;
	        writeedge(u,v,ofs);
            ofs<<endl;
        }
        ifs.close();
        ofs.close();
        
    }
     //delete temp file
    if(if_delete == 1){
        remove(readfilename.c_str());
    }

}



int main(int argc, char *argv[])
{

    if (argc != 6) {
        printf("input argc %d needed %d\n", argc, 6 );

        printf("usage: ./mergegraph infilename outfilename if_copy increase if_delete\n");
        return 0;
    }
    string ifile(argv[1]);
    string ofile(argv[2]);
    int if_copy = atoi(argv[3]);
    int increase = atoi(argv[4]);
    int if_delete = atoi(argv[5]);
    mergegraph(ifile, ofile ,if_copy, increase,if_delete);

    return 0;
}
