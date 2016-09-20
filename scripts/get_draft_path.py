#!/usr/bin/env python

import sys
import os
import subprocess
from parse_read import *
import numpy as np
import networkx as nx
import itertools
from pbcore.io import FastaIO




def rev_node(node):
    node_id = node.split('_')[0]
    return node_id + '_' + str(1-int(node.split('_')[1]))



def merge_nodes(g,in_node,out_node):

    weight = str(g.edge[in_node][out_node]['length'])

    if 'path' in g.node[in_node]:
        path1 = g.node[in_node]['path']
        weightspath1 = g.node[in_node]['weightspath']
    else:
        path1 = in_node
        weightspath1 = ''

    if 'path' in g.node[out_node]:
        path2 = g.node[out_node]['path']
        weightspath2 = ';' + g.node[out_node]['weightspath']
    else:
        path2 = out_node
        weightspath2 = ''

    g.node[in_node]['path'] = path1 + ';' + path2

    if weightspath1 == '':
        g.node[in_node]['weightspath'] = weight + weightspath2
    else:
        g.node[in_node]['weightspath'] = weightspath1 + ';' + weight + weightspath2

    for nodeB in g.successors(out_node):
        g.add_edge(in_node,nodeB,length=g.edge[out_node][nodeB]['length'])


    g.node[in_node]['cut_end'] = g.node[out_node]['cut_end']
    g.remove_node(out_node)



filedir = sys.argv[1]
filename = sys.argv[2]
graphml_path = sys.argv[3]

in_graph = nx.read_graphml(graphml_path)

# debug output
#for node in in_graph.nodes():
#    print node

#for edge in in_graph.edges():
#    print len(in_graph.edge[edge[0]][edge[1]])

reads = sorted(list(set([int(x.split("_")[0].lstrip("B")) for x in in_graph.nodes()])))

dbshow_reads = ' '.join([str(x+1) for x in reads])

DBshow_cmd = "DBshow "+ filedir+'/'+ filename+' '+dbshow_reads
stream = subprocess.Popen(DBshow_cmd.split(),
                                  stdout=subprocess.PIPE,bufsize=1)
reads_queried = parse_read(stream.stdout)
read_dict = {}
for read_id,read in itertools.izip(reads,reads_queried):
    rdlen = len(read[1])
#     print read
    read_dict[read_id] = read

# to simulate reads

# read_dict = {}
# for vertex in in_graph.nodes():
#     read_dict[int(vertex.split('_')[0])] = ['A','A'*50000]


complement = {'A':'T','C': 'G','T':'A', 'G':'C','a':'t','t':'a','c':'g','g':'c'}

# out_graphml_name = 'test.graphml'
out_graphml_name = filedir + '/' + filename +'_draft.graphml'

# outfile = 'test.edges.list'
outfile = filedir + '/' + filename + ".edges.list"


rev_comp_contig = True

out_graph = in_graph.copy()


# first we add some info to the graph for the cutting of contigs
for vert in out_graph.nodes():

    vert_id, vert_or = vert.split("_")
    vert_id = vert_id.lstrip("B")

    vert_len = len(read_dict[int(vert_id)][1])

    out_graph.node[vert]['cut_start'] = 0
    out_graph.node[vert]['cut_end'] = vert_len


    # SHOULD THIS USE THE RAW MATCHES?

    if out_graph.in_degree(vert) > 1:
        if vert_or == '0':
            out_graph.node[vert]['cut_start'] = max([out_graph.edge[x][vert]['read_b_match_start'] for x in out_graph.predecessors(vert)])
        else:
            out_graph.node[vert]['cut_start'] = vert_len - min([out_graph.edge[vert_id+'_0'][x]['read_a_match_start'] for x in out_graph.successors(vert_id+'_0')])

    if out_graph.out_degree(vert) > 1:

        if vert_or == '0':
            out_graph.node[vert]['cut_end'] = min([out_graph.edge[vert][x]['read_a_match_start'] for x in out_graph.successors(vert)])
        else:
            out_graph.node[vert]['cut_end'] = vert_len - max([out_graph.edge[x][vert_id+'_0']['read_b_match_start'] for x in out_graph.predecessors(vert_id+'_0')])




# next we merge the nodes in out_graph to form the contigs

nodes_to_merge = [x for x in out_graph.nodes() if out_graph.in_degree(x) == 1 and out_graph.out_degree(out_graph.predecessors(x)[0]) == 1]



# print len(read_dict[41260][1])
# print len(read_dict[4697][1])


while nodes_to_merge:

    cur_node = nodes_to_merge[0]

    prev_node = out_graph.predecessors(cur_node)[0]


    if prev_node != cur_node:
        merge_nodes(out_graph,prev_node,cur_node)
    else:
        out_graph.node[cur_node]['path'] = out_graph.node[cur_node]['path'] + ';' + cur_node
        out_graph.node[cur_node]['weightspath'] = out_graph.node[cur_node]['weightspath'] + ';' + str(out_graph.edge[prev_node][cur_node]['length'])
        out_graph.node[cur_node]['cut_end'] = len(read_dict[int(cur_node.split('_')[0])][1])


    nodes_to_merge.pop(0)

    # print len(nodes_to_merge)






# next we print the contigs out to the .edges.list file
contig_no = 0
print "Writing out_graph with "+str(len(out_graph.nodes()))+" contigs/nodes"


# we keep track of the already printed nodes so that reverse complement pairs are printed together
# we don't add to printed_nodes the "border" nodes so that we still have a partition of the nodes into contigs
printed_nodes = set()

# debug output

# for node in out_graph.nodes():
#     print node
#     print out_graph.node[node]


# for edge in out_graph.edges():
#     print edge
#     print out_graph.edge[edge[0]][edge[1]]




with open(outfile, 'w') as f:

    for vertex in out_graph.nodes():

        if rev_node(vertex) in printed_nodes:
            continue

        # single-node contig
        if 'path' not in out_graph.node[vertex]:
            f.write('>Unitig%d\n'%(contig_no))
            contig_no += 1

            # we repeat the same node twice so that the line is easily distinguishable (6 numbers)
            f.write('O %s %s %s %s %d %d\n'%(vertex.split('_')[0],vertex.split('_')[1]  , vertex.split('_')[0],
                vertex.split('_')[1], out_graph.node[vertex]['cut_start'], out_graph.node[vertex]['cut_end']) )

            printed_nodes = printed_nodes | set([vertex])

            f.write('>Unitig%d\n'%(contig_no))
            contig_no += 1

            vertex_rc = rev_node(vertex)
            f.write('O %s %s %s %s %d %d\n'%(vertex_rc.split('_')[0],vertex_rc.split('_')[1]  , vertex_rc.split('_')[0],
                vertex_rc.split('_')[1], out_graph.node[vertex_rc]['cut_start'], out_graph.node[vertex_rc]['cut_end']) )

            continue



        node_list = out_graph.node[vertex]['path'].split(';')
        weights_list = out_graph.node[vertex]['weightspath'].split(';')

        # print out_graph.node[vertex]['path']
        # print node_list
        # print out_graph.node[vertex]['weightspath']
        # print weights_list

        # print len(node_list),len(weights_list)

        if len(node_list) != len(weights_list)+1:
            print 'Something went wrong with contig '+str(contig_no)
            continue

        printed_nodes = printed_nodes | set(node_list)

        # print 'Unitig ' +str(contig_no)
        f.write('>Unitig%d\n'%(contig_no))
        contig_no += 1

        # prev_vert = out_graph.node[node_list[0]]['prev_node']
        # if prev_vert != '':

        if out_graph.in_degree(vertex) == 1 and out_graph.predecessors(vertex)[0] != vertex:

            prev_contig = out_graph.predecessors(vertex)[0]
            cut_start = out_graph.node[prev_contig]['cut_end']
            if out_graph.node[prev_contig].has_key('path'):
                nodeA = out_graph.node[prev_contig]['path'].split(';')[-1]
            else:
                nodeA = prev_contig

            nodeB = node_list[0]
            f.write('S %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0],
                nodeB.split('_')[1], out_graph.edge[prev_contig][vertex]['length'], cut_start) )
            nodeA = node_list[0]
            nodeB = node_list[1]
            f.write('T %s %s %s %s %s\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1], weights_list[0]) )

        else:
            nodeA = node_list[0]
            nodeB = node_list[1]
            f.write('S %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1],
                weights_list[0], out_graph.node[vertex]['cut_start']) )

        for i in range(1,len(weights_list)-1):
            nodeA = node_list[i]
            nodeB = node_list[i+1]
            f.write('T %s %s %s %s %s\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1], weights_list[i]) )

        if out_graph.out_degree(vertex) == 1 and out_graph.successors(vertex)[0] != vertex:

            nodeA = node_list[len(weights_list)-1]
            nodeB = node_list[len(weights_list)]
            f.write('T %s %s %s %s %s\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1], weights_list[-1]) )

            next_contig = out_graph.successors(vertex)[0]
            # we end this contig where the next one begins
            cut_end = out_graph.node[next_contig]['cut_start']

            nodeA = node_list[len(weights_list)]
            if out_graph.node[next_contig].has_key('path'):
                nodeB = out_graph.node[next_contig]['path'].split(';')[0]
            else:
                nodeB = next_contig

            f.write('E %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0],
                nodeB.split('_')[1], out_graph.edge[vertex][next_contig]['length'], cut_end) )

        else:

            nodeA = node_list[len(weights_list)-1]
            nodeB = node_list[len(weights_list)]
            f.write('E %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1],
                weights_list[-1], out_graph.node[vertex]['cut_end']) )




        # if we want reverse complement contigs, we print them next to each other

        if rev_comp_contig == False:
            continue


        f.write('>Unitig%d\n'%(contig_no))
        contig_no += 1



        if out_graph.out_degree(vertex) == 1 and out_graph.successors(vertex)[0] != vertex:

            next_contig = out_graph.successors(vertex)[0]

            nodeB = rev_node(node_list[len(weights_list)])
            if out_graph.node[next_contig].has_key('path'):
                nodeA = rev_node(out_graph.node[next_contig]['path'].split(';')[0])
            else:
                nodeA = rev_node(next_contig)

            # we start this contig where the previous (rc: next) one ended
            cut_start = len(read_dict[int(nodeA.split('_')[0])][1]) - out_graph.node[next_contig]['cut_start']

            f.write('S %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0],
                nodeB.split('_')[1], out_graph.edge[vertex][next_contig]['length'], cut_start) )

            nodeA = rev_node(node_list[len(weights_list)])
            nodeB = rev_node(node_list[len(weights_list)-1])
            f.write('T %s %s %s %s %s\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1], weights_list[-1]) )

        else:

            nodeA = rev_node(node_list[len(weights_list)])
            nodeB = rev_node(node_list[len(weights_list)-1])

            f.write('S %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1],
                weights_list[-1], len(read_dict[int(nodeA.split('_')[0])][1]) - out_graph.node[vertex]['cut_end']) )



        for i in range(len(weights_list)-1,1,-1):
            nodeA = rev_node(node_list[i])
            nodeB = rev_node(node_list[i-1])
            f.write('T %s %s %s %s %s\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1], weights_list[i-1]) )

        if out_graph.in_degree(vertex) == 1 and out_graph.predecessors(vertex)[0] != vertex:

            nodeA = rev_node(node_list[1])
            nodeB = rev_node(node_list[0])
            f.write('T %s %s %s %s %s\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1], weights_list[0]) )

            prev_contig = out_graph.predecessors(vertex)[0]

            nodeA = rev_node(node_list[0])

            if out_graph.node[prev_contig].has_key('path'):
                nodeB = rev_node(out_graph.node[prev_contig]['path'].split(';')[-1])
            else:
                nodeB = rev_node(prev_contig)
            cut_end = len(read_dict[int(nodeB.split('_')[0])][1]) - out_graph.node[prev_contig]['cut_end']

            f.write('E %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0],
                nodeB.split('_')[1], out_graph.edge[prev_contig][vertex]['length'], cut_end) )

        else:
            nodeB = rev_node(node_list[0])
            nodeA = rev_node(node_list[1])
            f.write('E %s %s %s %s %s %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], nodeB.split('_')[1],
                weights_list[0], len(read_dict[int(nodeB.split('_')[0])][1]) - out_graph.node[vertex]['cut_start']) )




print "Number of contigs: "+str(contig_no)

nx.write_graphml(out_graph,out_graphml_name)





