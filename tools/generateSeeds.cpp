// -*- coding: utf-8 -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
// generate Seed files
void generateSeeds(long long num, long long range){
    
    std::string output_name =std::string("seed_")+ std::to_string(num)+ std::string("_")+std::to_string(range)+std::string(".txt");
    std::ofstream ofs(output_name, std::ofstream::out);
    long long process = num;

    while (process>0) {
        std:: set<long long> setofseeds;
        long long u;
        u = rand()%range;
        auto it = setofseeds.find(u);
        if(it != setofseeds.end()){
            
        }else{
            setofseeds.insert(u);
            // read to file
            if( process != num){
                ofs << std::endl;
            }
             process--;
             ofs << u;
        }
       
    }
    printf("finished! \n");

}



int main(int argc, char *argv[])
{

    if (argc != 3) {
        printf("usage: ./generateSeeds seednum, range \n");
        return 0;
    }

    long num = atol(argv[1]);
    long range = atol(argv[2]);
    printf ("generate %d seeds, range 0 - %d \n",num, range);
    generateSeeds(num,range);


    return 0;
}
