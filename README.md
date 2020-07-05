# GLP
A light GPU system for variants of label propagation
Organization
--------

Code for "Graph Label Propagation on GPUs"

Compilation
--------

Requirements

* CMake &gt;= 2.8
* CUDA environment
* OpenMP

Compilation is done within the Root/ directory with CMake. 
Please make sure you have install CMake software compilation tool.
Configure cmakelist.txt appropriately before you start compile. 
To compile, simple rum

```
$ cmake .
$ make .
```

Running code in GLP
--------
Running code is also done within the Root/ directory. 
Use ".GPULP/ method iterations inputGraphPath outputPath ifDirectGraph" to run a label progagation for an input graph.
if "method" is set to 1, a classic label propagation will be applied to the input file. 
For layered Lp and SLPA, just set "method" to 2 and 3 respectively. If you want to implement your own version of LP.
Please overwrite APIs provide by us.

Examples
```
$ ./GPULP 1 20 ./datasets/text.bin output.bin 0
$ ./GPULP 2 20 ./datasets/text.bin output.bin 0
$ ./GPULP 3 20 ./datasets/text.bin output.bin 0
```
Input Format for GLP
--------

AdjacencyGraph  
0 1  
0 2  
1 3  
4 5  
6 7  

GLP support adacency graph of input. 
Each line represents one edge of the graph.
Note that input file should be in binary format.
Before you run LP algorithm, 
you can easily convert the input graph in "txt" format to "bin" format using "edgelist2bin" and "normalize" programs in the tools directory.

```
$ edgelist2bin test.txt
$ normalize test.bin
```

User-defined APIs
--------
**PickLabel(VertexId vid)**:Given a vertex vid, it decides vid’s label and write the label to the current label array L.  
**LoadNeighbor(VertexId vid)**:Given an edge vid, did, it returns the label as well as the label frequency for vertex did as a neighbor of vertex vid.  
**LabelScore(VertexId vid, LabelT l, double freq)**:Given a vertex $vid$, a label $l$ and  $l$'s frequency among $vid$'s neighbors, it returns a score of $l$ for $vid$.  
**UpdateVertex(VertexId vid, LabelT l, double score)**: Given a vertex $vid$, update the status of vertex $vid$ with label $l$ and $score$.  

you can overwrite those APIs, then run a customized LP after setting "method" to 4.
