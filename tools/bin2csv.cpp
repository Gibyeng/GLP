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
void bin2csv(std::string filename){
    long progress = 0;
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);
    ifs.seekg(0, std::ios::beg);
	const int bufsize = 2;
    int buf[bufsize];
    if(ifs.is_open()){
        std::cout << "find file..." <<std::endl;
    }
    std::string output_name = filename.substr(filename.find_last_of("/")+1,filename.find_last_of(".")-filename.find_last_of("/"))+std::string("csv");
    std::ofstream ofs(output_name, std::ofstream::out);
    while (ifs.read((char *) buf, sizeof(int) * bufsize)) {
        int u, v;
        u = buf[0];
		v = buf[1];
        if(progress!=0){
            ofs<<std::endl;
        }
        ofs << u <<','<< v;
        progress++;
        if ((progress & 0xffffff) == 0) {
            printf("%d \n", progress);
        }
    }
}



int main(int argc, char *argv[])
{

    if (argc != 2) {
        printf("usage: ./bin2csv filename\n");
        return 0;
    }

    const std::string filename(argv[1]);

    bin2csv(filename);


    return 0;
}
