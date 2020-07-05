
// -*- coding: utf-8 -*-
#pragma once

struct Test_result{
	//time	
	double load_graph_to_cpu = 0;
	double load_graph_to_gpu = 0;
	double other_proprocess = 0;
	double label_load = 0;
	double lp = 0;
	double tiny_node_process = 0;
	double small_node_process = 0;
	double medium_node_process = 0;
	double big_node_process = 0;
	double huge_node_process = 0;
	//graph distribution
	double tiny_distribution = 0;
	double small_distribution = 0;
	double medium_distribution = 0;
	double big_distribution = 0;
	double huge_distribution = 0;

};
