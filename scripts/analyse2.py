#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter



def de_clip(filename, n_iter):

    g = nx.MultiDiGraph()
    
    with open(filename,'r') as f:
        for line in f.xreadlines():
            l = line.strip().split()
            #print l2
            g.add_edge(l[0],l[1])
    
    print nx.info(g)
    degree_sequence=sorted(g.degree().values(),reverse=True)
    print Counter(degree_sequence)
    for i in range(n_iter):
        for node in g.nodes():
            if g.in_degree(node) == 0:
                g.remove_node(node)
    
        print nx.info(g)
    
    degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    print Counter(degree_sequence)
    
    
    
    
    
    try:
        import ujson
        mapping = ujson.load(open(filename.split('.')[0]+'.mapping.json'))
        
        print 'get mapping'
        
        for node in g.nodes():
            #print node
            if mapping.has_key(node):
                g.node[node]['aln_start'] = mapping[node][0]
                g.node[node]['aln_end'] = mapping[node][1]
                g.node[node]['aln_strand'] = mapping[node][2]
            else:
                g.node[node]['aln_start'] = 0
                g.node[node]['aln_end'] = 0
                g.node[node]['aln_strand'] = 0
                
    except:
        pass        
    
    nx.write_graphml(g, filename.split('.')[0]+'.graphml')
    
    print nx.number_weakly_connected_components(g)
    print nx.number_strongly_connected_components(g)
    
    
filename = sys.argv[1]
de_clip(filename, 20)
