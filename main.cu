#include <iostream>
#include <stdlib.h>

#include <cuda_profiler_api.h>

#include "header/myutil.h"
#include "header/graph.h"
#include "header/test_result.h"

#include "prop/layeredLprop.cuh"
#include "prop/SLPAprop.cuh"
#include "prop/Base5prop.cuh"

#include "prop/Test1prop.cuh"
#include "prop/classicLprop.cuh"
#include "prop/Test3prop.cuh"

int main(int argc, char *argv[])
 {

    //sorted graph using in method 15
    std::string graphpath {"./datasets/test.bin"};

	using Vertex = int32_t;
    using Edge = int32_t;

    //load graph
    using GraphT = CSRGraph<Vertex, Edge>;
    cudaSetDevice(0);
    cudaFree(0);
	// // switch method
	int method;
	int niter;
	std::string outputfile;
	std::string loadedgraphpath;

	CSRGraph<Vertex, Edge>* graph;
	if(argc == 1){
		method = 17;
		niter = 20;
		// graph, startvalue, if_direct
		// true - direct (14,15,16), false undirect(10-13,17,18)
		graph = new CSRGraph<Vertex, Edge>(graphpath,0,false);
	}else{
		if(argc == 6){
			method = atoi(argv[1]);
			niter = atoi(argv[2]);
			loadedgraphpath = argv[3];
			outputfile = argv[4];
			bool if_direct = (bool)atoi(argv[5]);
			graph = new CSRGraph<Vertex, Edge>(loadedgraphpath,0,if_direct);
		}else{
			printf("need 6 argc, input %d argc \n", argc);
			return 0;
		}
	}
	const Vertex n = graph->n;
	const Edge m = graph->m;
	std::cout <<"n = "<<n<<" m = "<<m<<std::endl;
	std::cout<<std::endl;
	double running_time = 0;
	Test_result* myresult = NULL;
    switch (method){

		case 1: //classic LP
		{
			auto propagator = new classicLprop<Vertex,Edge>(graph);
			myresult = propagator->run(niter);
			delete (propagator);
			break;
		}
		

		case 2:// layered LP  no optimization on loading neighbor label.
		{
			auto propagator = new MultiLprop<Vertex,Edge>(graph);
			running_time = propagator->run(niter);

			delete (propagator);
			break;
		}

		case 3:// SLPA
		{
			auto propagator = new SLPAprop<Vertex,Edge>(graph);
			running_time = propagator->run(niter);

			delete (propagator);
			break;
		}
		// user-defined lp
		case 4:{
			auto propagator = new UserdefinedLprop<Vertex,Edge>(graph);
			propagator->run(niter);

			delete (propagator);
			break;
		}


    }

    // free memory
    delete (graph);
    std::cout<<"success run";
    if(argc == 6){
    	std::ofstream ofs(outputfile, std::ofstream::out| std::ios::app);

    	if(myresult == NULL){
    		ofs <<"graph: "<< loadedgraphpath <<" method: "<< method << " iteration: "<< niter <<" time: "<< running_time ;
    		ofs<< std::endl;
    	}else{
    		ofs <<"graph: "<< loadedgraphpath <<" method: "<< method << " iteration: "<< niter <<" time: "<<" load graph to GPU: "<< myresult-> load_graph_to_gpu;
    		ofs<<" other proprocess: "<< myresult->other_proprocess << " label_load: " << myresult->label_load <<" LP: " <<  myresult->lp;
    		ofs<<" tiny_speed: "<< myresult->tiny_node_process << " small_speed: " << myresult->small_node_process <<" medium_speed: " <<myresult->medium_node_process << " big_speed: " <<myresult->big_node_process << " huge_speed: " << myresult-> huge_node_process;
    		ofs<<" tiny_distrbution: "<< myresult->tiny_distribution << " small_distribution: " << myresult->small_distribution <<" medium_distribution: " <<myresult->medium_distribution << " big_distribution: " <<myresult->big_distribution << " huge_distribution: " << myresult-> huge_distribution;
    		ofs<< std::endl;
    	}

    }

    return 0;
}
