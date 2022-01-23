import networkx as nx
import matplotlib.pyplot as plt
import sys
import random

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

true_pagerank = nx.pagerank(G)
true_pagerank = sorted(true_pagerank.items(), key=lambda x:x[1], reverse=True)
sample_pagerank = nx.pagerank(Sub_G)
#sample_pagerank = sorted(sample_pagerank.items(), key=lambda x:x[1], reverse=True)

sample_num = len(list(Sub_G.nodes))
sample_num = 10000
x = [i for i in range(sample_num)]
y = []
count = 0
top_nodes = []
top_k = []
for v, val in true_pagerank:
  top_nodes.append(v)
  y.append(float(val))
  if count < 1000:
    top_k.append(v)
  if count == sample_num-1:
    break
  count += 1

x_sample = []
y_sample = []
top_k_count = 0
count = 0
#for v, pr in sample_pagerank:
#  if v in top_nodes:
#    x_sample.append(count)
#    y_sample.append(float(pr))
#  if count < 1000:
#    if v in top_k:
#      top_k_count += 1
#  if count == sample_num-1:
#    break
#  count += 1
for v, pr in true_pagerank:
  if v in visited:
    x_sample.append(count)
    y_sample.append(float(sample_pagerank[v]))
  if count == sample_num-1:
    break
  count += 1

top_k_count = 0
for v in top_k:
  if v in visited:
    top_k_count += 1

print("top-k ratio " + str(top_k_count/1000))
fig, ax = plt.subplots()

normalized_y = []
for val in y:
  normalized_y.append(val/y[0])

normalized_sample_y = []
for val in y_sample:
  normalized_sample_y.append(val/y_sample[0])

ax.plot(x, normalized_y, color="blue")
ax.scatter(x_sample, normalized_sample_y, color="red", marker="x", s=10)
#ax.plot(x_sample, y_sample, color="red")
plt.show()
