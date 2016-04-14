#!/usr/bin/python

import networkx as nx
import random
import sys
from collections import Counter



# This script does a random condensation of the graph down to 2000 nodes

# python random_condensation.py ecoli.edges 2000

# It also keeps the ground truth on the graph through the condensation steps (if a json file is available)



def merge_path(g,in_node,node,out_node):
        
    g.add_edge(in_node,out_node,hinge_edge = -1,false_positive = 0)        
    g.remove_node(node)
    

def input1(flname):

    print "input1"

    g = nx.DiGraph()
    with open (flname) as f:
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

    print "input2"

    g = nx.DiGraph()
    with open (flname) as f:
        for lines in f:
            lines1=lines.split()
            #print lines1
            g.add_edge(lines1[0], lines1[1])   
    return g



def input3(flname):

    print "input3"
    # g = nx.DiGraph()
    g = nx.read_graphml(flname)


def de_clip(filename, n_nodes, hinge_list,gt_file):

    n_iter = 5

    
    f=open(filename)
    line1=f.readline()
    print line1
    f.close()

    extension = filename.split('.')[-1]

    if extension == 'graphml':
        g=input3(filename)
    elif len(line1.split()) !=2:
        g=input1(filename)
    else:
        g=input2(filename)

    
    print nx.info(g)
    degree_sequence=sorted(g.degree().values(),reverse=True)
    print Counter(degree_sequence)
    
    degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    print Counter(degree_sequence)
    
    try:
        import ujson
        mapping = ujson.load(open(gt_file))
        
        print 'getting mapping'
        mapped_nodes=0
        print str(len(mapping)) 
        print str(len(g.nodes()))
        for node in g.nodes():
            # print node
            node_base=node.split("_")[0]
            # print node_base

            #print node
            if mapping.has_key(node_base):
                g.node[node]['aln_start'] = min (mapping[node_base][0][0],mapping[node_base][0][1])
                g.node[node]['aln_end'] = max(mapping[node_base][0][1],mapping[node_base][0][0])
                g.node[node]['chr'] = mapping[node_base][0][2]
                mapped_nodes+=1
            else:
                # pass
                g.node[node]['aln_start'] = 0
                g.node[node]['aln_end'] = 0
                g.node[node]['aln_strand'] = 0


        for edge in g.edges_iter():
            in_node=edge[0]
            out_node=edge[1]
            # print  'akjdfakjhfakljh'
            if ((g.node[in_node]['aln_start'] < g.node[out_node]['aln_start'] and  
                g.node[out_node]['aln_start'] < g.node[in_node]['aln_end']) or 
                (g.node[in_node]['aln_start'] < g.node[out_node]['aln_end'] and
                g.node[out_node]['aln_end'] < g.node[in_node]['aln_end'])):
                g.edge[in_node][out_node]['false_positive']=0
            else:
                g.edge[in_node][out_node]['false_positive']=1

    except:
        raise
        # print "json "+filename.split('.')[0]+'.mapping.json'+" not found. exiting."
           
    print hinge_list

    print str(mapped_nodes)+" out of " +str(len(g.nodes()))+" nodes mapped."
    
    # for i in range(5):
    #     merge_simple_path(g)
    #     degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    #     print Counter(degree_sequence)

    in_hinges = set()
    out_hinges = set()
    num_iter=10000
    iter_done=0
    if hinge_list != None:
        print "Found hinge list."
        with open(hinge_list,'r') as f:
            for lines in f:
                lines1=lines.split()

                if lines1[2] == '1':
                  in_hinges.add(lines1[0]+'_0')
                  out_hinges.add(lines1[0]+'_1')
                elif lines1[2] == '-1':
                  in_hinges.add(lines1[0]+'_1')
                  out_hinges.add(lines1[0]+'_0')

        print str(len(in_hinges))+' hinges found.'

        for node in g.nodes():
            if node in in_hinges and node in out_hinges:
                g.node[node]['hinge']=100
            elif node in in_hinges:
                g.node[node]['hinge']=10
            elif node in out_hinges:
                g.node[node]['hinge']=-10
            else:
                g.node[node]['hinge']=0

        while len(g.nodes()) > n_nodes and iter_done < num_iter :
            node = g.nodes()[random.randrange(len(g.nodes()))]
            iter_done+=1
            # print iter_done
            if g.in_degree(node) == 1 and g.out_degree(node) == 1:

                base_node=node.split("_")[0]
                orintation = node.split("_")[1]
                # if orintation=='1':
                #     node2=base_node+'_0'
                # else:
                #     node2=base_node+'_1'

                # print node,node2

                in_node = g.in_edges(node)[0][0]
                out_node = g.out_edges(node)[0][1]

                if g.node[node]['hinge']==0 and g.node[in_node]['hinge']==0  and g.node[out_node]['hinge']==0:
                    if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                        if in_node != node and out_node != node and in_node != out_node:
                            bad_node=False
                            # print g.in_edges(node)
                            # print g.edge[g.in_edges(node)[0][0]][g.in_edges(node)[0][1]]
                            # print g.out_edges(node)
                            for in_edge in g.in_edges(node):
                                if g.edge[in_edge[0]][in_edge[1]]['false_positive']==1:
                                    bad_node=True
                            for out_edge in g.out_edges(node):
                                if g.edge[out_edge[0]][out_edge[1]]['false_positive']==1:
                                    bad_node=True 
                            if not bad_node:
                                #print in_node, node, out_node
                                merge_path(g,in_node,node,out_node)


                # print g.edge[edge1[0]][edge1[1]]['hinge_edge']

                for nd in g.nodes():
                    if len(nd.split("_"))==1:
                        print nd + " in trouble"
                # in_node = g.in_edges(node2)[0][0]
                # out_node = g.out_edges(node2)[0][1]
                # if g.node[node2]['hinge']==0 and g.node[in_node]['hinge']==0  and g.node[out_node]['hinge']==0:
                #     if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                #         if in_node != node2 and out_node != node2 and in_node != out_node:
                #             bad_node=False
                #             for in_edge in g.in_edges(node2):
                #                 if g.edge[in_edge]==1:
                #                     bad_node=True
                #             for out_edge in g.out_edges(node2):
                #                 if g.edge[out_edge]==1:
                #                     bad_node=True 
                #             if not bad_node:
                #                 #print in_node, node, out_node
                #                 merge_path(g,in_node,node2,out_node)


            # for nd in g.nodes():
            #     print nd

    else:
        while len(g.nodes()) > n_nodes:

            node = g.nodes()[random.randrange(len(g.nodes()))]



            if g.in_degree(node) == 1 and g.out_degree(node) == 1:

                # assert g.in_degree(node2) == 1 and g.out_degree(node2) == 1
                # edge_1 = g.out_edges(node)[0]
                # edge_2 = g.in_edges(node)[0]

                edge1 = g.out_edges(node)[0]
                edge2 = g.in_edges(node)[0]

                # print g.edge[edge1[0]][edge1[1]]['hinge_edge']

                if (g.edge[edge1[0]][edge1[1]]['hinge_edge'] == -1 and g.edge[edge2[0]][edge2[1]]['hinge_edge'] == -1):
                
                    in_node = g.in_edges(node)[0][0]
                    out_node = g.out_edges(node)[0][1]
                    if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                        if in_node != node and out_node != node and in_node != out_node:
                            #print in_node, node, out_node
                            merge_path(g,in_node,node,out_node)







    
    degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    print Counter(degree_sequence)

    
    nx.write_graphml(g, filename.split('.')[0]+'.sparse3.graphml')
    
    print nx.number_weakly_connected_components(g)
    print nx.number_strongly_connected_components(g)
  

if __name__ == "__main__":   
    filename = sys.argv[1]
    try :
        hinge_list=sys.argv[3]
        print "Found hinge list."
    except:
        hinge_list=None
        print "in except "+hinge_list

    de_clip(filename, int(sys.argv[2]),hinge_list, sys.argv[4])





