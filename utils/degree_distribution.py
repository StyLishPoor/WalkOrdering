import matplotlib.pyplot as plt
import networkx as nx
import sys

G = nx.read_edgelist(sys.argv[1], nodetype=int, comments='#')
#G = nx.gnp_random_graph(100, 0.5, directed=True)
m=3
#G = nx.barabasi_albert_graph(1000, m)

#degrees = [G.degree(n) for n in G.nodes()]
degree_freq = nx.degree_histogram(G)
degree_prob = list(map(lambda x:x/len(G.nodes()), degree_freq))
degrees = range(len(degree_freq))
plt.figure(figsize=(20, 20))
plt.loglog(degrees[m:], degree_prob[m:])
plt.show()
