// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
// read bin file then convert to csv file without any changes
// csv files for tigergraph
void bin2txt(std::string filename,bool if_direct){
    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    int buf[bufsize];
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
        if(!if_direct){
            ofs << v <<' '<< u;
        }
        progress++;
        if ((progress & 0xffffff) == 0) {
            printf("%d \n", progress);
        }
    }
}



int main(int argc, char *argv[])
{

    if (argc != 3) {
        printf("usage: ./bin2txt filename if_direct\n");
        return 0;
    }

    const std::string filename(argv[1]);
    const std::string directstr(argv[2]);
    bool if_direct = !(directstr=="0");
    if(if_direct){
        std::cout<< "direct: keep edges unchange"<<std::endl;

    }else{
        std::cout<< "undirect: double edges, only used in ligra"<<std::endl;
    }
    bin2txt(filename,if_direct);


    return 0;
}
