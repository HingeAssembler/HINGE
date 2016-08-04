#!/usr/bin/env python

import networkx as nx
import sys
from collections import Counter

def longest_path(G):
    dist = {} # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0]+1,v) for v in G.pred[node]] 
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node,(length,_)  = max(dist.items(), key=lambda x:x[1])
    path = []
    while length > 0:
        path.append(node)
        length,node = dist[node]
    return list(reversed(path))
    
    
filename = sys.argv[1]


g = nx.DiGraph()

with open(filename,'r') as f:
    for line in f.xreadlines():
        g.add_edge(*(line.strip().split('->')))


print nx.info(g)
degree_sequence=sorted(nx.degree(g).values(),reverse=True)
print Counter(degree_sequence)

for i in range(7):
    for node in g.nodes():
        if g.in_degree(node) == 0:
            g.remove_node(node)
            
    print nx.info(g)

degree_sequence=sorted(nx.degree(g).values(),reverse=True)
print Counter(degree_sequence)


def rev(string):
    if string[-1] == '\'':
        return string[:-1]
    else:
        return string+'\''

for edge in g.edges():
    g.add_edge(rev(edge[1]), rev(edge[0]))
    #print edge
    #print rev(edge[1]), rev(edge[0])

print nx.info(g)
nx.write_graphml(g, filename.split('.')[0]+'.graphml')
#print(list(nx.dfs_edges(g,sys.argv[2])))
#p=nx.shortest_path(g)
