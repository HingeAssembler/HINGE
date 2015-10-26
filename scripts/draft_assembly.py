#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter

def linearize(filename):
    graph_name = filename.split('.')[0]+'.graphml'
    g = nx.read_graphml(graph_name)
    
    print nx.info(g)
    
    # get first strong connected component
    
    con = list(nx.strongly_connected_component_subgraphs(g))
    
    con.sort(key = lambda x:len(x), reverse = True)
    print [len(item) for item in con]
    
    print nx.info(con[0])
    
    with open(filename.split('.')[0]+'.linear.edges', 'w') as f:
        for item in nx.dfs_preorder_nodes(con[0]):
            f.write(item + '\n')

filename = sys.argv[1]
linearize(filename)