#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter

def merge_simple_path(g):
    for node in g.nodes():
        #print g.in_degree(node), g.out_degree(node)
        if g.in_degree(node) == 1 and g.out_degree(node) == 1:
            
            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    #print in_node, node, out_node
                    merge_path(g,in_node,node,out_node)
                    
                              
def merge_two_nodes(g):
    for node in g.nodes():
        if g.in_degree(node) == 1 and g.out_degree(node) == 0:
            in_node = g.in_edges(node)[0][0]
            if g.out_degree(in_node) == 1:
                if in_node != node:
                    node_id = g.graph['aval']
                    g.graph['aval'] += 1
                    g.add_node(str(node_id), 
                        count = g.node[in_node]['count'] + g.node[node]['count'],
                        read = g.node[in_node]['read'] + '_' + g.node[node]['read'],
                        #aln_chr = g.node[node]['aln_chr']
                        )
                    g.remove_node(in_node)
                    g.remove_node(node)
            

def merge_path(g,in_node,node,out_node):
    #ov1 = find_overlap(g.node[in_node]['bases'], g.node[node]['bases'])
    #ov2 = find_overlap(g.node[node]['bases'], g.node[out_node]['bases'])
    
    node_id = g.graph['aval']
    g.graph['aval'] += 1
    #length = g.node[node]['length'] + g.node[in_node]['length'] + g.node[out_node]['length'] - ov1 - ov2
    #cov = (g.node[in_node]['cov'] * g.node[in_node]['length'] + g.node[node]['cov'] * g.node[node]['length']  + \
    #g.node[out_node]['cov'] * g.node[out_node]['length'])/float(length)
    #bases = g.node[in_node]['bases'][:-ov1] + g.node[node]['bases'] + g.node[out_node]['bases'][ov2:]
    
    g.add_node(str(node_id), 
        count = g.node[in_node]['count'] + g.node[node]['count'] + g.node[out_node]['count'],
        read = g.node[in_node]['read'] + '_' +  g.node[node]['read'] + '_' +g.node[out_node]['read'],
        #aln_chr = g.node[node]['aln_chr']
    )
    #g.add_node(str(node_id)+'-', bases = reverse_comp_bases(bases), length = length, cov = cov)
    
    for edge in g.in_edges(in_node):
        g.add_edge(edge[0],str(node_id))
    
    for edge in g.out_edges(out_node):
        g.add_edge(str(node_id),edge[1])
    
        
    g.remove_node(in_node)
    g.remove_node(node)
    g.remove_node(out_node)
    
def input1(flname):
    g = nx.DiGraph()
    with open (filename) as f:
        for lines in f:
            lines1=lines.split()
            #print lines1
            if len(lines1) < 5:
                continue
            #print lines1
            g.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4], hinge_edge=int(lines1[5]))
            g.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))
    return g
            
def input2(flname):
    g = nx.DiGraph()
    with open (filename) as f:
        for lines in f:
            lines1=lines.split()
            #print lines1
            g.add_edge(lines1[0], lines1[1])   
    return g

def run(filename, n_iter):
    
    
    f=open(filename)
    line1=f.readline()
    print line1
    f.close()
    if len(line1.split()) !=2:
	g=input1(filename)
    else:
	g=input2(filename)
    
    
    
    print nx.info(g)
    
    
    for node in g.nodes():
        g.node[node]['count'] = 1
        g.node[node]['read'] = node
        
        
    degree_sequence=sorted(g.degree().values(),reverse=True)
    print Counter(degree_sequence)
    for i in range(n_iter):
        for node in g.nodes():
            if g.in_degree(node) == 0:
                g.remove_node(node)
    
        print nx.info(g)
        degree_sequence=sorted(nx.degree(g).values(),reverse=True)
        print Counter(degree_sequence)
    
    degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    print Counter(degree_sequence)
    
    
    g.graph['aval'] = 1000000000
    
    for i in range(5):
        merge_simple_path(g)
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
    
    nx.write_graphml(g, filename.split('.')[0]+'_condensed.graphml')
    
    print nx.number_weakly_connected_components(g)
    print nx.number_strongly_connected_components(g)
    
    
filename = sys.argv[1]
run(filename, 5)
