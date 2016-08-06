#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter


# This script condenses the graph down, creates a gfa with for the condensed graph, and computes the contig N50

# python condense_graph_create_gfa_compute_n50.py ecoli.edges

# The conditions in lines 23 and 24 are meant to prevent nodes corresponding to different strands to be merged 
# (and should be commented out if this is not desired, or if a json is not available)


def merge_simple_path(g):
    for node in g.nodes():
        if g.in_degree(node) == 1 and g.out_degree(node) == 1:
            
            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node: # and in_node != out_node:
                    # if g.node[in_node]['aln_strand']==g.node[node]['aln_strand'] or max(g.node[in_node]['aln_strand'],g.node[node]['aln_strand']) == 5:
                        # if g.node[out_node]['aln_strand']==g.node[node]['aln_strand'] or max(g.node[out_node]['aln_strand'],g.node[node]['aln_strand']) == 5:
                    #print in_node, node, out_node
                    merge_path(g,in_node,node,out_node)


# def merge_path(g,in_node,node,out_node):

#     if g.edge[in_node][node]['intersection'] == 1 and g.edge[node][out_node]['intersection'] == 1:
#         g.add_edge(in_node,out_node,hinge_edge = -1,intersection = 1,z=0)
#     else:
#         g.add_edge(in_node,out_node,hinge_edge = -1,intersection = 0,z=0)

#     g.remove_node(node)



def merge_path_old(g,in_node,node,out_node):
    #ov1 = find_overlap(g.node[in_node]['bases'], g.node[node]['bases'])
    #ov2 = find_overlap(g.node[node]['bases'], g.node[out_node]['bases'])
    
    node_id = g.graph['aval']
    g.graph['aval'] += 1
    #length = g.node[node]['length'] + g.node[in_node]['length'] + g.node[out_node]['length'] - ov1 - ov2
    #cov = (g.node[in_node]['cov'] * g.node[in_node]['length'] + g.node[node]['cov'] * g.node[node]['length']  + \
    #g.node[out_node]['cov'] * g.node[out_node]['length'])/float(length)
    #bases = g.node[in_node]['bases'][:-ov1] + g.node[node]['bases'] + g.node[out_node]['bases'][ov2:]

    overlap1 = g.edge[in_node][node][0]['overlap']
    overlap2 = g.edge[node][out_node][0]['overlap']

    length0 = g.node[in_node]['length']
    length1 = g.node[node]['length']
    length2 = g.node[out_node]['length']


    if overlap1 > min(length0,length1):
        print "problem here:"
        print overlap1, length0, length1


    g.add_node(str(node_id),length = length0+length1+length2 - overlap1 - overlap2)
    #g.add_node(str(node_id)+'-', bases = reverse_comp_bases(bases), length = length, cov = cov)
    
    for cur_edge in g.in_edges(in_node):

        # print g.edge[cur_edge[0]][cur_edge[1]][0]['overlap']

        g.add_edge(cur_edge[0],str(node_id),overlap = g.edge[cur_edge[0]][cur_edge[1]][0]['overlap'])
    
    for cur_edge in g.out_edges(out_node):
        g.add_edge(str(node_id),cur_edge[1],overlap = g.edge[cur_edge[0]][cur_edge[1]][0]['overlap'])
        
    g.remove_node(in_node)
    g.remove_node(node)
    g.remove_node(out_node)
    



def merge_path(g,in_node,node,out_node):
    #ov1 = find_overlap(g.node[in_node]['bases'], g.node[node]['bases'])
    #ov2 = find_overlap(g.node[node]['bases'], g.node[out_node]['bases'])
    
    node_id = g.graph['aval']
    g.graph['aval'] += 1
    #length = g.node[node]['length'] + g.node[in_node]['length'] + g.node[out_node]['length'] - ov1 - ov2
    #cov = (g.node[in_node]['cov'] * g.node[in_node]['length'] + g.node[node]['cov'] * g.node[node]['length']  + \
    #g.node[out_node]['cov'] * g.node[out_node]['length'])/float(length)
    #bases = g.node[in_node]['bases'][:-ov1] + g.node[node]['bases'] + g.node[out_node]['bases'][ov2:]

    overlap1 = g.edge[in_node][node][0]['overlap']
    overlap2 = g.edge[node][out_node][0]['overlap']

    length0 = g.node[in_node]['length']
    length1 = g.node[node]['length']
    length2 = g.node[out_node]['length']


    if overlap1 > min(length0,length1):
        print "problem here:"
        print overlap1, length0, length1

    g.node[in_node]['length'] = length0 + length1 - overlap1

    # g.add_edge(in_node,out_node,overlap = g.edge[cur_edge[0]][cur_edge[1]][0]['overlap'])
    g.add_edge(in_node,out_node,overlap = g.edge[node][out_node][0]['overlap'])        

    # g.remove_node(in_node)
    g.remove_node(node)
    # g.remove_node(out_node)



def comp_n50(contig_vec):
    if len(contig_vec) == 0:
        return 0
    sorted_lengths = sorted(contig_vec)

    total_length = sum(contig_vec)

    half_length = 0.5*total_length

    min_n50 = sorted_lengths[-1]
    max_n50 = 0

    for i in range(len(sorted_lengths)):
    #if len(sorted_lengths) % 2 == 0:
    #  sum_1 = sum(sorted_lengths[0:i])
    #  sum_2 = sum(sorted_lengths[i:])
    #else:
    #  sum_1 = sum(sorted_lengths[0:i+1])
    #  sum_2 = sum(sorted_lengths[i:])
        sum_1 = sum(sorted_lengths[0:i+1])
        sum_2 = sum(sorted_lengths[i:])
        if sum_1 >= half_length and sum_2 >= half_length:
            min_n50 = min(sorted_lengths[i],min_n50)
            max_n50 = max(sorted_lengths[i],max_n50)

  # print "Min N50: "+str(min_n50)
  # print "Max N50: "+str(max_n50)

    return 0.5*(min_n50+max_n50)



def de_clip(filename, n_iter):

    print "version 0.4"

    g = nx.MultiDiGraph()
    
    # count = 0

    with open(filename,'r') as f:
        for line in f.xreadlines():
            l = line.strip().split()
            #print l2
            g.add_edge(l[0],l[1],overlap=int(l[2])/2)
            # if count < 10:
            #     print l[0], l[1], l[2]
            #     count += 1

            # node0start = int(l[7][1:])
            # node0end = int(l[8][:-1])

            # g.node[l[0]]['length'] = node0end - node0start
            g.node[l[0]]['length'] = int(l[3])

            # node1start = int(l[9][1:])
            # node1end = int(l[10][:-1])

            # g.node[l[1]]['length'] = node1end - node1start
            g.node[l[1]]['length'] = int(l[4])
    
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
    
    
    g.graph['aval'] = 1000000000
    
    for i in range(8):
        merge_simple_path(g)
        degree_sequence=sorted(nx.degree(g).values(),reverse=True)
        print Counter(degree_sequence)
    
       
    
    nx.write_graphml(g, filename.split('.')[0]+'.graphml')
    
    print nx.number_weakly_connected_components(g)
    print nx.number_strongly_connected_components(g)


    # Next we create the gfa file


    outputfile = filename.split('.')[0]+'.gfa'
    with open(outputfile, 'w') as fout:

        for cur_node in g.nodes():

            node_length = g.node[cur_node]['length']
            node_str = 'A'*node_length
            node_str = node_str + '\n'

            fout.write("NODE "+str(cur_node)+' 0 0 0 0 0\n')
            fout.write(node_str)
            fout.write(node_str)
            # print "NODE "+str(node)

        for arc in g.edges():
            fout.write("ARC "+str(arc[0])+' '+str(arc[1])+' 0\n')



    # Compute N50

    contig_lengths = []

    for cur_node in g.nodes():
        contig_lengths.append(g.node[cur_node]['length'])

    print sorted(contig_lengths)

    print "N50 = "+str(comp_n50(contig_lengths))

    print [x for x in contig_lengths if x > 40000]

    print sum([x for x in contig_lengths if x <= 40000])

    # print "N50 = "+str(comp_n50([x for x in contig_lengths if x > 100000]))

    print "number of contigs (> 40000): "+str(len([x for x in contig_lengths if x > 40000]))
    
    
filename = sys.argv[1]
de_clip(filename, 5)
