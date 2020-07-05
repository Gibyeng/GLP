// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>

void convert_tsv(std::string filename){
    long progress = 0;
    std::ifstream ifs(filename, std::ios::ate);
    ifs.seekg(0, std::ios::beg);
    std::string line;
    if(ifs.is_open()){
        std::cout << "find file..." <<std::endl;
    }
    int c_num;
    int i_num;
    std:: vector<long long > c_list;
    std:: vector<long long > i_list;
    // first read...
    std::cout << "first read begin..." <<std::endl;
    std::set<long long> c_set;
    std::set<long long> i_set;
    while (getline(ifs,line)) {
        std::stringstream ss(line);
        std::string nodestr;
        long long u, v;
        getline(ss,nodestr,'\t');
        u = stol(nodestr);
        getline(ss,nodestr,'\t');
        v = stol(nodestr);
        c_set.insert(u);
        i_set.insert(v);
        progress++;
         if ((progress & 0xffffff) == 0) {
            printf("%d \n", progress);
        }
        
    }
    progress = 0;
    c_list.assign (c_set.begin(),c_set.end());
    i_list.assign (i_set.begin(),i_set.end());
    c_num = c_list.size();
    i_num = i_list.size(); 
    std::cout << "number of u: "<<c_num<<std::endl;
    std::cout << "number of v: "<<i_num<<std::endl;
    // second read...
    std::cout << "second read begin..." <<std::endl;
    std::string output_name_txt = std::string("conv_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("txt");
    std::string output_name_bin = std::string("conv_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");
    std::ofstream ofs_txt(output_name_txt, std::ofstream::out);
    std::ofstream ofs_bin(output_name_bin, std::ofstream::out);
    int edge_num = 0;
    long long old_u = 0;
    long long old_v = 0;
    ifs.close();
    ifs.open(filename, std::ios::ate);
    ifs.seekg(0, std::ios::beg);
    while (getline(ifs,line)) {
        std::stringstream ss(line);
        std::string nodestr;
        long long u, v;
        getline(ss,nodestr,'\t');
        u = stol(nodestr);
        getline(ss,nodestr,'\t');
        v = stol(nodestr);
        // remove dup edges
        if(u==old_u && v == old_v){
            continue;
        }
        edge_num++;
        old_u = u;
        old_v = v;
        //find u 
        auto it  = std::find(c_list.begin(),c_list.end(),u);
        int u_r = std::distance(c_list.begin(), it);
        //find v
        it  = std::find(i_list.begin(),i_list.end(), v);
        int v_r = std::distance(i_list.begin(), it);
        v_r += c_num;
        if(progress!=0){
            ofs_txt<<std::endl;
            ofs_bin<<std::endl;
        }
        ofs_txt << u_r <<' '<< v_r;
        ofs_bin << u_r <<' '<< v_r;
        progress++;
        if ((progress & 0xffffff) == 0) {
            printf("%d \n", progress);
        }
    }
    printf("edge number: %d \n",edge_num);
}



int main(int argc, char *argv[])
{

    if (argc != 2) {
        printf("usage: ./convert_tsv filename\n");
        return 0;
    }

    const std::string filename(argv[1]);

    convert_tsv(filename);


    return 0;
}
