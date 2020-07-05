#include <iostream>
#include <stdlib.h>

#include <cuda_profiler_api.h>

#include "myutil.h"
#include "graph.h"

#include "prop/Lprop.cuh"
#include "prop/CMLprop.cuh"
#include "prop/MultiLprop.cuh"
#include "prop/BasicLprop.cuh"
#include "prop/Base1prop.cuh"
#include "prop/Base2prop.cuh"
#include "prop/Base3prop.cuh"
#include "prop/PushLprop.cuh"
#include "prop/MultiGprop.cuh"
#include "prop/MBase1prop.cuh"
#include "prop/MBase2prop.cuh"
#include "prop/CPUlp.cuh"
#include "prop/seedLprop.cuh"
#include "prop/layeredLprop.cuh"
#include "prop/Base4prop.cuh"
#include "prop/SLPAprop.cuh"
#include "prop/Base5prop.cuh"
#include "test_result.h"
#include "prop/Test1prop.cuh"
#include "prop/Test2prop.cuh"
#include "prop/Test3prop.cuh"

int main(int argc, char *argv[])
 {
//	// start with 1
//    std::string graphpath0 {"/home/changye/Documents/GPULP/datasets/sorted_direct_test.bin"};
////    std::string graphpath0 {"/home/changye/Documents/GPULP/datasets/test.bin"};
//    std::string graphpath1 {"/home/changye/Documents/GPULP/datasets/com-amazon.ungraph.bin"};
//    std::string graphpath2 {"/home/changye/Documents/GPULP/datasets/com-dblp.ungraph.bin"};
//   // std::string graphpath20 {"/home/changye/Documents/GPULP/datasets/normalized_com-amazon.ungraph.txt"};
//    std::string graphpath3 {"/home/changye/Documents/GPULP/datasets/com-youtube.ungraph.bin"};
//    //std::string graphpath4 {"/home/changye/Documents/GPULP/datasets/com-friendster.ungraph.bin"};
//    // start with 0
//    std::string graphpath4 {"/home/changye/Documents/GPULP/datasets/roadNet-CA.bin"};
//    std::string graphpath5 {"/home/changye/Documents/GPULP/datasets/com-lj.ungraph.bin"};
//    std::string graphpath6 {"/home/changye/Documents/GPULP/datasets/com-orkut.ungraph.bin"};
//    std::string graphpath7 {"/home/changye/Documents/GPULP/datasets/wiki-en.bin"};
//    std::string graphpath8 {"/home/changye/Documents/GPULP/datasets/com-friendster.ungraph.txt.bin"};
//    std::string graphpath9 {"/home/changye/Documents/GPULP/datasets/uk-2002.bin"};

    //sorted graph using in method 15
    std::string graphpath10 {"/home/changye/Documents/GPULP/datasets/sorted_com-youtube.ungraph.bin"};
    std::string graphpath11 {"/home/changye/Documents/GPULP/datasets/sorted_com-lj.ungraph.bin"};
    std::string graphpath12 {"/home/changye/Documents/GPULP/datasets/sorted_com-dblp.ungraph.bin"};
    std::string graphpath13 {"/home/changye/Documents/GPULP/datasets/sorted_com-orkut.ungraph.bin"};
    std::string graphpath14 {"/home/changye/Documents/GPULP/datasets/sorted_direct_wikipedia_en.bin"};
    std::string graphpath15 {"/home/changye/Documents/GPULP/datasets/sorted_direct_uk-2002.bin"};
    std::string graphpath16 {"/home/changye/Documents/GPULP/datasets/sorted_direct2_twitter_rv.bin"};
    std::string graphpath18 {"/home/changye/Documents/GPULP/datasets/sorted_roadNet-CA.bin"};
    std::string graphpath17 {"/home/changye/Documents/GPULP/datasets/sorted_conv_aligraph.bin"};

    //man made database
    std::string graphpath19{"/home/changye/Documents/GPULP/datasets/sorted_graph500_19_4.bin"};
	std::string graphpath20{"/home/changye/Documents/GPULP/datasets/sorted_graph500_20_4.bin"};
	std::string graphpath21{"/home/changye/Documents/GPULP/datasets/sorted_graph500_21_4.bin"};
	std::string graphpath22{"/home/changye/Documents/GPULP/datasets/sorted_graph500_22_4.bin"};
	std::string graphpath23{"/home/changye/Documents/GPULP/datasets/sorted_graph500_23_4.bin"};

    std::string graphpath29{"/home/changye/Documents/GPULP/datasets/sorted_graph500_19_8.bin"};
    std::string graphpath30{"/home/changye/Documents/GPULP/datasets/sorted_graph500_20_8.bin"};
    std::string graphpath31{"/home/changye/Documents/GPULP/datasets/sorted_graph500_21_8.bin"};
    std::string graphpath32{"/home/changye/Documents/GPULP/datasets/sorted_graph500_22_8.bin"};
    std::string graphpath33{"/home/changye/Documents/GPULP/datasets/sorted_graph500_23_8.bin"};
    std::string graphpath34{"/home/changye/Documents/GPULP/datasets/sorted_graph500_24_8.bin"};
    std::string graphpath35{"/home/changye/Documents/GPULP/datasets/sorted_graph500_25_8.bin"};

    std::string graphpath39{"/home/changye/Documents/GPULP/datasets/sorted_graph500_19_16.bin"};
	std::string graphpath40{"/home/changye/Documents/GPULP/datasets/sorted_graph500_20_16.bin"};
	std::string graphpath41{"/home/changye/Documents/GPULP/datasets/sorted_graph500_21_16.bin"};
	std::string graphpath42{"/home/changye/Documents/GPULP/datasets/sorted_graph500_22_16.bin"};
	std::string graphpath43{"/home/changye/Documents/GPULP/datasets/sorted_graph500_23_16.bin"};
    std::string graphpath44{"/home/changye/Documents/GPULP/datasets/sorted_graph500_24_16.bin"};
	std::string graphpath45{"/home/changye/Documents/GPULP/datasets/sorted_graph500_25_16.bin"};

	std::string graphpath46{"/home/changye/Documents/GPULP/datasets/sorted_graph500_23_12.bin"};
	std::string graphpath47{"/home/changye/Documents/GPULP/datasets/sorted_graph500_23_20.bin"};
	std::string graphpath48{"/home/changye/Documents/GPULP/datasets/sorted_graph500_23_24.bin"};

	std::string graphpath50{"/home/changye/Documents/GPULP/datasets/D/sorted_D2_part1.bin"};
	std::string graphpath51{"/home/changye/Documents/GPULP/datasets/D/D7_part1.bin"};
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
		graph = new CSRGraph<Vertex, Edge>(graphpath10,0,false);
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

        case 1:// label propagation final
        {
        	auto propagator = new Lprop<Vertex,Edge>(graph);
        	running_time = propagator->run(niter);

			delete (propagator);
			break;
        }
        case 2:// warp + block (two kernel) baseline 1
		{
			auto propagator = new Base1prop<Vertex,Edge>(graph);
			propagator->run(niter);

			delete (propagator);
			break;
		}
        case 3:// warp + block + shared memory baseline 2
        {
        	auto propagator = new Base2prop<Vertex,Edge>(graph);
			propagator->run(niter);

			delete (propagator);
			break;
        }

        case 4://block baseline 3 -- hash table (Ghash)
	   {
			auto propagator = new Base3prop<Vertex,Edge>(graph);
			propagator->run(niter);

			delete (propagator);
			break;
	   }
        case 5: // push-base method update
        {
			auto propagator = new PushLprop<Vertex,Edge>(graph);
			propagator->run(niter);

			delete (propagator);
			break;
	    }
        case 6: // multiGPU
		{
			auto propagator = new MultiGprop<Vertex,Edge>(graph);
			running_time = propagator->run(niter);

			delete (propagator);
			break;
		}

        case 7:// warp + block (two kernel) baseline 1 multiGPU
		{
			auto propagator = new MBase1prop<Vertex,Edge>(graph);
			propagator->run(niter);

			delete (propagator);
			break;
		}

	case 8:// warp + block + shared memory baseline 2 multiGPU
	{
		auto propagator = new MBase2prop<Vertex,Edge>(graph);
		propagator->run(niter);

		delete (propagator);
		break;
	}

	case 9:{
		auto propagator = new CPUlp<Vertex,Edge>(graph);
//		running_time = propagator->sync_run_1(niter);
//		running_time = propagator->layered_run(niter);
		running_time = propagator->slpa_run(niter);
		delete (propagator);
		break;
	}

	// seed label prop	
	case 10:{
		auto propagator = new seedLprop<Vertex,Edge>(graph);
		propagator->run(niter);

		delete (propagator);
		break;
	}


	case 11:// layered LP  no optimization on loading neighbor label.
	{
		auto propagator = new MultiLprop<Vertex,Edge>(graph);
		running_time = propagator->run(niter);

		delete (propagator);
		break;
	}
	case 12:// layered LP GPU hashtable
	{
		auto propagator = new Base4prop<Vertex,Edge>(graph);
		running_time = propagator->run(niter);

		delete (propagator);
		break;
	}
	case 13:// SLPA
	{
		auto propagator = new SLPAprop<Vertex,Edge>(graph);
		running_time = propagator->run(niter);

		delete (propagator);
		break;
	}
	case 14:// SLPA hash_table
	{
		auto propagator = new Base5prop<Vertex,Edge>(graph);
		running_time = propagator->run(niter);

		delete (propagator);
		break;
	}

	case 15: // performance test ~ Lprop- notuse
	{
		auto propagator = new Test1prop<Vertex,Edge>(graph);
		myresult = propagator->run(niter);
		delete (propagator);
		break;
	}
	case 16: // performance test ~ Lprop without label loading opt
	{
		auto propagator = new Test2prop<Vertex,Edge>(graph);
		myresult = propagator->run(niter);
		delete (propagator);
		break;
	}

	case 17: // performance test ~ Lprop without label loading opt
	{
		auto propagator = new Test3prop<Vertex,Edge>(graph);
		myresult = propagator->run(niter);
		delete (propagator);
		break;
	}
//        case 3: // layered LP: apply optimization 1
//        {
//        	auto propagator = new MultiOptLprop<Vertex,Edge>(graph);
//			propagator->run(niter);
//			delete (propagator);
//			break;
//        }

//        case 4:// cm-based
//	   {
//			auto propagator = new CMLprop<Vertex,Edge>(graph);
//			propagator->run(niter);
//
//			delete (propagator);
//			break;
//	   }

//		case 5:// hash base line
//	   {
//			auto propagator = new BasicLprop<Vertex,Edge>(graph);
//			propagator->run(niter);
//
//			delete (propagator);
//			break;
//	   }
//		case 6: // experiment test
//		{
//			auto propagator = new MultiOptLprop<Vertex,Edge>(graph);
//			//(int niter,bool opt_load,bool onestep,int changeiter,int delta,int alpha,std::string filename)
//			propagator->experiment(199,true,true,10,1000,10,std::string{"1.csv"});
//			delete (propagator);
//			break;
//		}

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
