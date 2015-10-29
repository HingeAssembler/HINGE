#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter

def de_clip(filename, n_iter):

    g = nx.DiGraph()
    
    with open(filename,'r') as f:
        for line in f.xreadlines():
            l = line.strip().split(',')
            l2 = l[0].split('->') + [l[1]]
            print l2
            g.add_edge(l2[0],l2[1],weight = int(l2[2])/100000.0, ew = int(l2[2]))
    
    
    print nx.info(g)
    degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    print Counter(degree_sequence)
    
    for i in range(n_iter):
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
    print nx.info(g)
    nx.write_graphml(g, filename.split('.')[0]+'.graphml')
    
    
filename = sys.argv[1]
de_clip(filename, 30)
    