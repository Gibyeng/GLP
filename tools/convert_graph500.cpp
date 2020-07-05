// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>

void convert_graph500(std::string filename){
    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    long long buf[bufsize];
    if(ifs.is_open()){
        std::cout << "find file..." <<std::endl;
    }
    std::string output_name_txt = std::string("conv_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("txt");
    std::string output_name_bin = std::string("conv_")+filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("bin");
    std::ofstream ofs_txt(output_name_txt, std::ofstream::out);
    std::ofstream ofs_bin(output_name_bin, std::ofstream::out);
    while (ifs.read((char *) buf, sizeof(long long) * bufsize)) {
        int u, v;
        u = abs((int)buf[0]);
		v = abs((int)buf[1]);
        if(progress!=0){
            ofs_txt<<std::endl;
            ofs_bin<<std::endl;
        }
        ofs_txt << u <<' '<< v;
        ofs_bin << u <<' '<< v;
        progress++;
        if ((progress & 0xffffff) == 0) {
            printf("%d \n", progress);
        }
    }
}

void convert_bin_edgelist_32(std::string filename){
    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    long int buf[bufsize];
    if(ifs.is_open()){
        std::cout << "find file..." <<std::endl;
    }
    std::string output_name = filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("txt");
    std::ofstream ofs(output_name, std::ofstream::out);
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        int u, v;
        u = buf[0];
		v = buf[1];
        if(progress!=0){
            ofs<<std::endl;
        }
        ofs << u <<' '<< v;
        progress++;
        if ((progress & 0xffffff) == 0) {
            printf("%d \n", progress);
        }
    }
}

int main(int argc, char *argv[])
{

    if (argc != 2) {
        printf("usage: ./bin2edgelist filename\n");
        return 0;
    }

    const std::string filename(argv[1]);

    convert_graph500(filename);


    return 0;
}
