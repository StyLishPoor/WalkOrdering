import sys
import networkx as nx

G = nx.read_edgelist(sys.argv[1], comments='#')
G = G.to_directed()

nx.write_edgelist(G, sys.argv[2], data=False)
