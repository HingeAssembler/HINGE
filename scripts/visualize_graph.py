#!/usr/bin/python

import networkx as nx
import sys

filename = sys.argv[1]


g = nx.DiGraph()

with open(filename,'r') as f:
    for line in f.xreadlines():
        g.add_edge(*(line.strip().split('->')))


print nx.info(g)
nx.write_graphml(g, filename.split('.')[0]+'.graphml')
#print(list(nx.dfs_edges(g,sys.argv[2])))
