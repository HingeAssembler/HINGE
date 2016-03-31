#!/usr/bin/python

import networkx as nx
import random
import sys
from collections import Counter



# This script does a random condensation of the graph down to 2000 nodes

# python random_condensation.py ecoli.edges 2000

# It also keeps the ground truth on the graph through the condensation steps (if a json file is available)



def merge_simple_path(g):
    for node in g.nodes():
        if g.in_degree(node) == 1 and g.out_degree(node) == 1:
            
            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    #print in_node, node, out_node
                    merge_path(g,in_node,node,out_node)

def merge_path(g,in_node,node,out_node):
    #ov1 = find_overlap(g.node[in_node]['bases'], g.node[node]['bases'])
    #ov2 = find_overlap(g.node[node]['bases'], g.node[out_node]['bases'])
    
    node_id = g.graph['aval']
    g.graph['aval'] += 1
    #length = g.node[node]['length'] + g.node[in_node]['length'] + g.node[out_node]['length'] - ov1 - ov2
    #cov = (g.node[in_node]['cov'] * g.node[in_node]['length'] + g.node[node]['cov'] * g.node[node]['length']  + \
    #g.node[out_node]['cov'] * g.node[out_node]['length'])/float(length)
    #bases = g.node[in_node]['bases'][:-ov1] + g.node[node]['bases'] + g.node[out_node]['bases'][ov2:]
    g.add_node(str(node_id),aln_end = g.node[out_node]['aln_end'])
    #g.add_node(str(node_id)+'-', bases = reverse_comp_bases(bases), length = length, cov = cov)
    
    for edge in g.in_edges(in_node):
        g.add_edge(edge[0],str(node_id))
    
    for edge in g.out_edges(out_node):
        g.add_edge(str(node_id),edge[1])
    
        
    g.remove_node(in_node)
    g.remove_node(node)
    g.remove_node(out_node)
    
    

def de_clip(filename, n_nodes):

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
            if g.degree(node) < 2:
                g.remove_node(node)
    
        print nx.info(g)
        degree_sequence=sorted(nx.degree(g).values(),reverse=True)
        print Counter(degree_sequence)
    
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

    
    g.graph['aval'] = 1000000000
    
    # for i in range(5):
    #     merge_simple_path(g)
    #     degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    #     print Counter(degree_sequence)



    while len(g.nodes()) > n_nodes:
        node = g.nodes()[random.randrange(len(g.nodes()))]

        if g.in_degree(node) == 1 and g.out_degree(node) == 1:
            
            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    #print in_node, node, out_node
                    merge_path(g,in_node,node,out_node)



    
     
    
    nx.write_graphml(g, filename.split('.')[0]+'.sparse.graphml')
    
    print nx.number_weakly_connected_components(g)
    print nx.number_strongly_connected_components(g)
    
    
filename = sys.argv[1]
de_clip(filename, sys.argv[2])

