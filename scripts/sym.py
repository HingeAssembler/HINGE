#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter
def rev(s):
    if s[-1] == '\'':
        return s[:-1]
    else:
        return s+'\''

def merge(filename):
    with open(filename+'.1','r') as f:
        lines1 = [' '.join(item.strip().split()[:2]) for item in f.readlines()]
    with open(filename+'.2','r') as f:
        lines2 = [' '.join(item.strip().split()[:2]) for item in f.readlines()]
        
    #print lines1
    #print len(lines1)
    #print len(lines2)
    
    
    #for item in (set(lines1) - set(lines2)):
    #    print item

    #print len(set(lines1) - set(lines2))
    
    g1 = nx.DiGraph()
    for line in lines1:
        l = line.strip().split()
        g1.add_edge(l[0],l[1])
    
    #for i in range(30):
    #    for node in g1.nodes():
    #        if g1.in_degree(node) == 0:
    #            g1.remove_node(node)
    
    
    g = nx.DiGraph()
    
    for line in set(lines1) & set(lines2):
        l = line.strip().split()
        #print l2
        g.add_edge(l[0],l[1])
    
    for edge in g1.edges():
        if g.has_node(edge[0]) and g.out_degree(edge[0]) != 0:
            pass
        else:
            g.add_edge(edge[0],edge[1])
            g.add_edge(rev(edge[1]),rev(edge[0]))
            
    #for i in range(30):
    #    for node in g.nodes():
    #        if g.in_degree(node) == 0:
    #            g.remove_node(node)
            
            
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
merge(filename)