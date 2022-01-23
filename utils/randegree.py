import networkx as nx
import matplotlib.pyplot as plt
import sys
import random
import numpy as np

G = nx.read_edgelist(sys.argv[1], comments='#', nodetype=int)
G = G.to_directed()

max_degree = -1
max_id = -1

for v in list(G.nodes):
  if G.degree[v] > max_degree:
    max_degree = G.degree[v]
    max_id = v

start_node = max_id
current_node = start_node
visited = []
visited.append(start_node)

while (len(visited) < len(G.nodes) * float(sys.argv[2])):
  next_node = random.choice(list(G.neighbors(current_node)))
  if (random.random() <= float(G.degree[current_node])/G.degree[next_node]):
    if next_node not in visited:
      visited.append(next_node)
    if G.degree[next_node] == 0:
      current_node = start_node
    else:
      current_node = next_node
  else:
    continue
    
Sub_G = G.subgraph(visited)

for i in range(1, int(sys.argv[3])):
  print("-------[" + str(i) + "]-------")
  sub_list = np.array_split(visited, i)

  for group in sub_list:
    average_degree = 0;
    for v in group:
      average_degree += Sub_G.degree[v]
    print(average_degree/len(Sub_G.nodes))
