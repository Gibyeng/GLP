import networkit as nk
import os
print("reading graph")
G = nk.graphio.EdgeListReader(' ',0,commentPrefix='%',continuous=False,directed=False).read("./tools/conv_graph500_23_24.txt")
n = G.numberOfNodes()
m = G.numberOfEdges()
print(n, m)
LP = nk.community.PLP(G).run()
print(LP.getTiming())
print(LP.numberOfIterations())
print(sum(LP.getTiming())/1000)