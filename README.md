# GLP
A light GPU system for variants of label propagation

Organization
--------
Code for "Graph Label Propagation on GPUs"

Download
--------
There is a public tool to download the code from anonymous4open. (https://github.com/ShoufaChen/clone-anonymous4open)

Compilation
--------

Requirements

* CMake &gt;= 2.8
* CUDA environment
* OpenMP

Compilation is done within the Root/ directory with CMake. 
Please make sure you have installed CMake software compilation tool.
Configure cmakelist.txt appropriately before you start compile. 
To compile, simple run:

```
$ cmake .
$ make .
```

Running code in GLP
--------
Running code is also done within the Root/ directory. 
Use ".GPULP/ method iterations inputGraphPath outputPath ifDirectedGraph ifWeightedGraph" to run a label progagation for an input graph.
if "method" is set to 1, a classic label propagation will be applied to the input file. 
For layered Lp and SLPA, just set "method" to 2 and 3 respectively. If you want to implement your own version of LP.
Please overwrite APIs provide by us.

Examples
```
$ ./GPULP 1 20 ./datasets/text.bin output.bin 0 0
or
$ ./GPULP method=1 iter=20 input=./datasets/text.bin output=output.bin ifdirect=0 ifweighted=0
```
Input Format for GLP
--------

AdjacencyGraph(unweighted)  
0 1  
0 2  
1 3  
4 5  
6 7  
AdjacencyGraph(weighted)  
v 0 1  
v 1 2  
v 2 3  
e 0 1 1  
e 1 2 1  

GLP support adjacency graph of input. For unweighted graph,
each line represents one edge of the graph.
For weighted graph, we use "v vertexid vertexLabel" and "e vertexid vertexid edgeLabel" format.
Note that input file should be in binary format.
Before you run LP algorithm, 
you can easily convert the input graph in "txt" format to "bin" format using "edgelist2bin" and "normalize" programs in the tools/ directory.

```
$ edgelist2bin test.txt
$ normalize test.bin
```
User-defined APIs
--------
**PickLabel(VertexId vid)**:Given a vertex vid, it decides vid’s label and write the label to the current label array L.  
**LoadNeighbor(VertexId vid)**:Given an edge vid, did, it returns the label as well as the label frequency for vertex did as a neighbor of vertex vid.  
**LabelScore(VertexId vid, LabelT l, double freq)**:Given a vertex vid, a label $l$ and  $l$'s frequency among vid's neighbors, it returns a score of l for vid.  
**UpdateVertex(VertexId vid, LabelT l, double score)**: Given a vertex vid, update the status of vertex vid with label l and score.  

you can overwrite those APIs in file Root/header/kernel.cuh, then run a customized LP after setting "method" to 4.

Example of classic LP
```c++
__device__ void PickLabel(VertexIdvid){
    /*copy Lnextto L */
    G-> Attr-> L[vid] = G-> Attr-> Lnext[vid];
}

__device__ pair<LabelT, double> LoadNeighbor(VertexIdvid, VertexIdsid){
    LabelTl = G-> Attr-> L[sid];
    /*accumulates frequency of l by one */
    return pair<double,LabelT>(l, 1.0);
}

__device__ double LabelScore(VertexIdvid, LabelTl, double freq){
    /*return freqas its label score*/
    return freq;
}

__device__ void UpdateVertex(VertexIdvid, LabelTl, double freq){
    return;
}
```

You can find more implements of various LP in the tutorial/ directory.