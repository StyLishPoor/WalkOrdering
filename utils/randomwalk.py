import sys
import networkx as nx
import random

G = nx.read_edgelist(sys.argv[1], comments='#', nodetype=int)
G = G.to_directed()

print(G.degree, key=lambda x: x[1], reverse=True)

#print(nx.number_of_nodes(G))
#graph_nodes = list(G.nodes)
#graph_nodes.sort()
#maxid=graph_nodes[-1]
#idlist = [i for i in range(maxid+1)]
#random.shuffle(idlist)
#
#
#output_graph = sys.argv[2]
#with open(output_graph, mode='w') as f:
#  for edge in G.edges:
#    #print(edge, maxid)
#    s = str(idlist[edge[0]]) + " " + str(idlist[edge[1]]) + '\n'
#    f.write(s)
#
