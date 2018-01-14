#!/usr/bin/env python

# coding: utf-8

# In[115]:

import networkx as nx
import random
import sys
import os
import numpy as np
import ujson
from colormap import rgb2hex
import operator
import matplotlib.colors
import configparser

import datetime


# In[3]:



def write_graph(G,flname):
    with open(flname,'w') as f:
        for edge in G.edges_iter():
            f.write(str(edge[0])+'\t'+str(edge[1])+'\n')


# In[4]:

def write_graph2(G,Ginfo,flname):

    count_no = 0
    count_yes = 0

    with open(flname,'w') as f:
        for edge in G.edges_iter():


            if (edge[0],edge[1]) not in Ginfo:
                count_no += 1
                print "not found"
                continue
            else:
                count_yes += 1


#             line = Ginfo[(edge[0],edge[1])]
#             line_sp = line.split(' ')

#             f.write(str(edge[0])+' '+str(edge[1]))
#             for j in range(2,len(line_sp)):
#                 f.write(' '+line_sp[j])

            f.write(Ginfo[(edge[0],edge[1])]+'\n')

    print count_no, count_yes







# In[7]:

def prune_graph(graph,in_hinges,out_hinges,reverse=False):

    H=nx.DiGraph()
    if reverse:
        G=nx.reverse(graph,copy=True)
    else:
        G=graph
    start_nodes = [x for x in G.nodes() if G.in_degree(x) ==0]

    in_hinges = list(in_hinges.intersection(set(G.nodes())))
    out_hinges = list(out_hinges.intersection(set(G.nodes())))

    if reverse:
        for node in in_hinges:
            for successor in G.successors(node):
#                 H.add_edge(node,successor)
                H.add_node(successor)
        for node in out_hinges:
            H.add_node(node)
    else:
        for node in out_hinges:
            for successor in G.successors(node):
#                 H.add_edge(node,successor)
                H.add_node(successor)
        for node in in_hinges:
            H.add_node(node)
    map(H.add_node,start_nodes)
    all_vertices=set(G.nodes())
    current_vertices=set(H.nodes())
    undiscovered_vertices=all_vertices-current_vertices
    last_discovered_vertices=current_vertices
    while undiscovered_vertices:
        discovered_vertices_set=set([x for node in last_discovered_vertices
                                  for x in G.successors(node)
                                  if x not in current_vertices])
        for vertex in discovered_vertices_set:
            for v_predecessor in G.predecessors(vertex):
                if v_predecessor in current_vertices:
                    H.add_edge(v_predecessor,vertex)
                    break
        current_vertices=current_vertices.union(discovered_vertices_set)
#         print len(undiscovered_vertices)
        if len(discovered_vertices_set)==0:
            print last_discovered_vertices
            print 'did not reach all nodes'
            print 'size of G: '+str(len(G.nodes()))
            print 'size of H: '+str(len(H.nodes()))
#             return H

            rand_node = list(undiscovered_vertices)[0]

            discovered_vertices_set.add(rand_node)


        last_discovered_vertices=discovered_vertices_set
        undiscovered_vertices=all_vertices-current_vertices
#     if reverse:
#         for vertex in out_hinges:
#             for v_predecessor in G.predecessors(vertex):
#                 H.add_edge(v_predecessor,vertex)
#     else:
#         for vertex in in_hinges:
#             for v_predecessor in G.predecessors(vertex):
#                 H.add_edge(v_predecessor,vertex)
    if reverse:
        for node in in_hinges:
            for successor in G.successors(node):
                H.add_edge(node,successor)
        for node in out_hinges:
            for predecessor in G.predecessors(node):
                H.add_edge(predecessor,node)
    else:
        for node in out_hinges:
            for successor in G.successors(node):
                H.add_edge(node,successor)
        for node in in_hinges:
            for predecessor in G.predecessors(node):
                H.add_edge(predecessor,node)
    if reverse:
        return nx.reverse(H)
    return H


# In[8]:

def dead_end_clipping(G,threshold):
#     H=nx.DiGraph()
    H = G.copy()
    start_nodes = set([x for x in H.nodes() if H.in_degree(x) ==0])

    for st_node in start_nodes:
        cur_path = [st_node]

        if len(H.successors(st_node)) == 1:
            cur_node = H.successors(st_node)[0]
            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1 and len(cur_path) < threshold + 2:
                cur_path.append(cur_node)
                cur_node = H.successors(cur_node)[0]


        if len(cur_path) <= threshold:
            for vertex in cur_path:
                H.remove_node(vertex)

    end_nodes = set([x for x in H.nodes() if H.out_degree(x) ==0])

    for end_node in end_nodes:
        cur_path = [end_node]
        if len(H.predecessors(end_node)) == 1:
            cur_node = H.predecessors(end_node)[0]
            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1 and len(cur_path) < threshold + 2:
                cur_path.append(cur_node)
                cur_node = H.predecessors(cur_node)[0]

        if len(cur_path) <= threshold:
            for vertex in cur_path:
                H.remove_node(vertex)

    return H



def rev_node(node):
    node_id = node.split('_')[0]

    return node_id + '_' + str(1-int(node.split('_')[1]))


def dead_end_clipping_sym(G,threshold,print_debug = False):

    H = G.copy()
    
#     print("Searching for dead-ends")
    start_nodes = set([x for x in H.nodes() if H.in_degree(x) ==0])

    num_start_nodes = len(start_nodes)
    st_node_count = 0
    
    for st_node in start_nodes:
        
#         if st_node_count % clipping_pc == 0:
#             update_progress(st_node_count/clipping_pc)        

        st_node_count += 1
        sys.stdout.write('\r[clip] Processing dead-ends: {}/{}'.format(st_node_count,num_start_nodes))
        sys.stdout.flush()
    
        if H.has_node(st_node) == False:
            continue

        cur_path = [st_node]

        cur_node = st_node
        if print_debug:
            print '----0'
            print st_node

        if len(H.successors(st_node)) == 1:
            cur_node = H.successors(st_node)[0]
            
            if print_debug:
                print '----1'

            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1 and len(cur_path) < threshold + 2:
                cur_path.append(cur_node)

                if print_debug:
                    print cur_node

                cur_node = H.successors(cur_node)[0]

                if len(cur_path) > threshold + 1:
                    break


        if print_debug:
            print '----2'
            print cur_path

        if len(cur_path) <= threshold and (H.in_degree(cur_node) > 1 or H.out_degree(cur_node) == 0):
            for vertex in cur_path:
                # try:
                if print_debug:
                    print 'about to delete ',vertex,rev_node(vertex)
                H.remove_node(vertex)
                H.remove_node(rev_node(vertex))
                # except:
                    # pass
                if print_debug:
                    print 'deleted ',vertex,rev_node(vertex)


    print ""
    return H



# In[9]:


# This function is no longer used. See z_clipping_sym
def z_clipping(G,threshold,in_hinges,out_hinges,print_z = False):
    H = G.copy()

    start_nodes = set([x for x in H.nodes() if H.out_degree(x) > 1 and x not in out_hinges])

    for st_node in start_nodes:
        for sec_node in H.successors(st_node):

            if H.out_degree(st_node) == 1:
                break

            cur_node = sec_node
            cur_path = [[st_node,cur_node]]

            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:
                cur_path.append([cur_node,H.successors(cur_node)[0]])
                cur_node = H.successors(cur_node)[0]

                if len(cur_path) > threshold + 1:
                    break

            if len(cur_path) <= threshold and H.in_degree(cur_node) > 1 and H.out_degree(st_node) > 1 and cur_node not in in_hinges:
                if print_z:
                    print cur_path

                for edge in cur_path:
                    H.remove_edge(edge[0],edge[1])
                for j in range(len(cur_path)-1):
                    H.remove_node(cur_path[j][1])

    end_nodes = set([x for x in H.nodes() if H.in_degree(x) > 1 and x not in in_hinges])

    for end_node in end_nodes:
        for sec_node in H.predecessors(end_node):

            if H.in_degree(end_node) == 1:
                break


            cur_node = sec_node
            cur_path = [[cur_node,end_node]]

            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:
                cur_path.append([H.predecessors(cur_node)[0],cur_node])
                cur_node = H.predecessors(cur_node)[0]

                if len(cur_path) > threshold + 1:
                    break

            if len(cur_path) <= threshold and H.out_degree(cur_node) > 1 and H.in_degree(end_node) > 1 and cur_node not in out_hinges:
                if print_z:
                    print cur_path
                for edge in cur_path:
                    H.remove_edge(edge[0],edge[1])
                for j in range(len(cur_path)-1):
                    H.remove_node(cur_path[j][0])

    return H



def z_clipping_sym(G,threshold,in_hinges,out_hinges,print_z = False):

    H = G.copy()
    G0 = G.copy()

    start_nodes = set([x for x in H.nodes() if H.out_degree(x) > 1 and x not in out_hinges])

    for st_node in start_nodes:

        try:  # need this because we are deleting nodes inside loop
            H.successors(st_node)
        except:
            continue

        for sec_node in H.successors(st_node):

            if H.out_degree(st_node) == 1:
                break

            cur_node = sec_node
            cur_path = [[st_node,cur_node]]

            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:

                cur_path.append([cur_node,H.successors(cur_node)[0]])
                cur_node = H.successors(cur_node)[0]

                if len(cur_path) > threshold + 1:
                    break

            if len(cur_path) <= threshold and H.in_degree(cur_node) > 1 and H.out_degree(st_node) > 1 and cur_node not in in_hinges:
                if print_z:
                    print cur_path

                for edge in cur_path:

                    G0.edge[edge[0]][edge[1]]['z'] = 1
                    G0.edge[rev_node(edge[1])][rev_node(edge[0])]['z'] = 1

                    try:
                        H.remove_edge(edge[0],edge[1])
                        H.remove_edge(rev_node(edge[1]),rev_node(edge[0]))

                    except:
                        pass

                for j in range(len(cur_path)-1):

                    G0.node[cur_path[j][1]]['z'] = 1
                    G0.node[rev_node(cur_path[j][1])]['z'] = 1

                    try:
                        H.remove_node(cur_path[j][1])
                        H.remove_node(rev_node(cur_path[j][1]))

                    except:
                        pass


    return H, G0






# In[48]:

def merge_path(g,in_node,node,out_node):

    g.add_edge(in_node,out_node,hinge_edge = -1,false_positive = 0, z=0)

    g.remove_node(node)



# In[121]:

def random_condensation(G,n_nodes,check_gt = False):

    g = G.copy()

    max_iter = 20000
    iter_cnt = 0

    while len(g.nodes()) > n_nodes and iter_cnt < max_iter:

        iter_cnt += 1

        node = g.nodes()[random.randrange(len(g.nodes()))]

        if g.in_degree(node) == 1 and g.out_degree(node) == 1:

            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    #print in_node, node, out_node
#                     merge_path(g,in_node,node,out_node)

                    bad_node=False
                    if check_gt:
                        for in_edge in g.in_edges(node):
                            if g.edge[in_edge[0]][in_edge[1]]['false_positive']==1:
                                bad_node=True
                        for out_edge in g.out_edges(node):
                            if g.edge[out_edge[0]][out_edge[1]]['false_positive']==1:
                                bad_node=True
                    if not bad_node:
                        #print in_node, node, out_node
                        merge_path(g,in_node,node,out_node)

    if iter_cnt >= max_iter:
        print "couldn't finish sparsification"+str(len(g.nodes()))

    return g



def random_condensation_sym(G,n_nodes,check_gt = False):

    g = G.copy()

    max_iter = 20000
    iter_cnt = 0

    while len(g) > n_nodes and iter_cnt < max_iter:

        iter_cnt += 1

        node = g.nodes()[random.randrange(len(g))]

        if g.in_degree(node) == 1 and g.out_degree(node) == 1:

            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    #print in_node, node, out_node
#                     merge_path(g,in_node,node,out_node)

                    bad_node=False
                    if check_gt:
                        for in_edge in g.in_edges(node):
                            if g.edge[in_edge[0]][in_edge[1]]['false_positive']==1:
                                bad_node=True
                        for out_edge in g.out_edges(node):
                            if g.edge[out_edge[0]][out_edge[1]]['false_positive']==1:
                                bad_node=True
                    if not bad_node:
                        #print in_node, node, out_node

                        try:
                            merge_path(g,in_node,node,out_node)
                            merge_path(g,rev_node(out_node),rev_node(node),rev_node(in_node))
                        except:
                            pass

    if iter_cnt >= max_iter:
        print "[clip] Sparsification ended with {} nodes.".format(str(len(g)))

    return g


# In[118]:

def random_condensation2(g,n_nodes):

    g = G.copy()


    max_iter = 20000
    iter_cnt = 0

    while len(g.nodes()) > n_nodes and iter_cnt < max_iter:

        iter_cnt += 1

        node = g.nodes()[random.randrange(len(g.nodes()))]

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



    if iter_cnt >= max_iter:
        print "couldn't finish sparsification: "+str(len(g.nodes()))


    return g




def bubble_bursting_sym(H,threshold,print_bubble = False):

    start_nodes = set([x for x in H.nodes() if H.out_degree(x) == 2])

    for st_node in start_nodes:

        try:  # need this because we are deleting nodes inside loop
            H.successors(st_node)[1]
        except:
            continue

        sec_node = H.successors(st_node)[0]

        cur_node = sec_node
        cur_path = [[st_node,cur_node]]

        while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:

            cur_path.append([cur_node,H.successors(cur_node)[0]])
            cur_node = H.successors(cur_node)[0]

            if len(cur_path) > threshold + 1:
                break

        end_node0 = cur_node
        cur_node = H.successors(st_node)[1]
        alt_path = [[st_node,cur_node]]

        while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:

            alt_path.append([cur_node,H.successors(cur_node)[0]])
            cur_node = H.successors(cur_node)[0]

            if len(alt_path) > threshold + 1:
                break


        if len(cur_path) <= threshold and len(alt_path) <= threshold and end_node0 == cur_node:

            if print_bubble:
                print 'found bubble'

            for edge in cur_path:

                # try:
                H.remove_edge(edge[0],edge[1])
                H.remove_edge(rev_node(edge[1]),rev_node(edge[0]))

                # except:
                #     pass

            for j in range(len(cur_path)-1):

                # try:
                H.remove_node(cur_path[j][1])
                H.remove_node(rev_node(cur_path[j][1]))

                # except:
                #     pass


    return H


def resolve_rep(g,rep_path,in_node,out_node):

    prefix = 'B'

    g.add_edge(in_node,prefix + rep_path[0],
        length=g.edge[in_node][rep_path[0]]['length'],
        read_a_match_start=g.edge[in_node][rep_path[0]]['read_a_match_start'],
        read_a_match_end=g.edge[in_node][rep_path[0]]['read_a_match_end'],
        read_b_match_start=g.edge[in_node][rep_path[0]]['read_b_match_start'],
        read_b_match_end=g.edge[in_node][rep_path[0]]['read_b_match_end'],
        read_a_match_start_raw=g.edge[in_node][rep_path[0]]['read_a_match_start_raw'],
        read_a_match_end_raw=g.edge[in_node][rep_path[0]]['read_a_match_end_raw'],
        read_b_match_start_raw=g.edge[in_node][rep_path[0]]['read_b_match_start_raw'],
        read_b_match_end_raw=g.edge[in_node][rep_path[0]]['read_b_match_end_raw'])
    g.remove_edge(in_node,rep_path[0])

    g.add_edge(prefix+rep_path[-1],out_node,
               length=g.edge[rep_path[-1]][out_node]['length'],
        read_a_match_start=g.edge[rep_path[-1]][out_node]['read_a_match_start'],
        read_a_match_end=g.edge[rep_path[-1]][out_node]['read_a_match_end'],
        read_b_match_start=g.edge[rep_path[-1]][out_node]['read_b_match_start'],
        read_b_match_end=g.edge[rep_path[-1]][out_node]['read_b_match_end'],
        read_a_match_start_raw=g.edge[rep_path[-1]][out_node]['read_a_match_start_raw'],
        read_a_match_end_raw=g.edge[rep_path[-1]][out_node]['read_a_match_end_raw'],
        read_b_match_start_raw=g.edge[rep_path[-1]][out_node]['read_b_match_start_raw'],
        read_b_match_end_raw=g.edge[rep_path[-1]][out_node]['read_b_match_end_raw'])
    g.remove_edge(rep_path[-1],out_node)


    g.add_edge(rev_node(prefix + rep_path[0]),rev_node(in_node),
               length =g.edge[rev_node(rep_path[0])][rev_node(in_node)]['length'],
        read_a_match_start=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_a_match_start'],
        read_a_match_end=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_a_match_end'],
        read_b_match_start=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_b_match_start'],
        read_b_match_end=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_b_match_end'],
        read_a_match_start_raw=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_a_match_start_raw'],
        read_a_match_end_raw=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_a_match_end_raw'],
        read_b_match_start_raw=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_b_match_start_raw'],
        read_b_match_end_raw=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_b_match_end_raw'])
    g.remove_edge(rev_node(rep_path[0]),rev_node(in_node))
    g.add_edge(rev_node(out_node),rev_node(prefix+rep_path[-1]),
               length=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['length'],
        read_a_match_start=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_a_match_start'],
        read_a_match_end=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_a_match_end'],
        read_b_match_start=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_b_match_start'],
        read_b_match_end=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_b_match_end'],
        read_a_match_start_raw=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_a_match_start_raw'],
        read_a_match_end_raw=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_a_match_end_raw'],
        read_b_match_start_raw=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_b_match_start_raw'],
        read_b_match_end_raw=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_b_match_end_raw'])
    g.remove_edge(rev_node(out_node),rev_node(rep_path[-1]))




    for i in range(0,len(rep_path)-1):
        g.add_edge(prefix+rep_path[i],prefix+rep_path[i+1],
                   length=g.edge[rep_path[i]][rep_path[i+1]]['length'],
            read_a_match_start=g.edge[rep_path[i]][rep_path[i+1]]['read_a_match_start'],
            read_a_match_end=g.edge[rep_path[i]][rep_path[i+1]]['read_a_match_end'],
            read_b_match_start=g.edge[rep_path[i]][rep_path[i+1]]['read_b_match_start'],
            read_b_match_end=g.edge[rep_path[i]][rep_path[i+1]]['read_b_match_end'],
            read_a_match_start_raw=g.edge[rep_path[i]][rep_path[i+1]]['read_a_match_start_raw'],
            read_a_match_end_raw=g.edge[rep_path[i]][rep_path[i+1]]['read_a_match_end_raw'],
            read_b_match_start_raw=g.edge[rep_path[i]][rep_path[i+1]]['read_b_match_start_raw'],
            read_b_match_end_raw=g.edge[rep_path[i]][rep_path[i+1]]['read_b_match_end_raw'])
        g.add_edge(rev_node(prefix+rep_path[i+1]),rev_node(prefix+rep_path[i]),
                   length =g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['length'],
            read_a_match_start=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_a_match_start'],
            read_a_match_end=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_a_match_end'],
            read_b_match_start=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_b_match_start'],
            read_b_match_end=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_b_match_end'],
            read_a_match_start_raw=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_a_match_start_raw'],
            read_a_match_end_raw=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_a_match_end_raw'],
            read_b_match_start_raw=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_b_match_start_raw'],
            read_b_match_end_raw=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_b_match_end_raw'])




def loop_resolution(g,max_nodes,flank,print_debug = False):

    starting_nodes =  [x for x in g.nodes() if g.out_degree(x) == 2]

    if print_debug:
        print '----'
        print starting_nodes

    tandem = []

    for st_node in starting_nodes:


        if g.out_degree(st_node) != 2:
            continue

        if print_debug:
            print '----'
            print st_node

        loop_len = 0

        for first_node in g.successors(st_node):


            if g.out_degree(st_node) != 2:
                continue

            if print_debug:
                print '----'
                print first_node

            other_successor = [x for x in g.successors(st_node) if x != first_node][0]

            next_node = first_node
            if print_debug:
                print 'going on loop'

            loop_len = 0
            prev_edge = g[st_node][next_node]
            node_cnt = 0
            while g.in_degree(next_node) == 1 and g.out_degree(next_node) == 1 and node_cnt < max_nodes:
                node_cnt += 1
                in_node = next_node
                next_node = g.successors(next_node)[0]
                loop_len += abs(g[in_node][next_node]['read_a_match_start'] - prev_edge['read_b_match_start'])
                prev_edge = g[in_node][next_node]

            if node_cnt >= max_nodes:
                continue

            if print_debug:
                print "length in loop " + str(loop_len)
            len_in_loop = loop_len
            first_node_of_repeat = next_node

            if g.in_degree(next_node) == 2:
                prev_node = [x for x in g.predecessors(next_node) if x != in_node][0]
                node_cnt = 0
                while g.in_degree(prev_node) == 1 and g.out_degree(prev_node) == 1:
                    node_cnt += 1
                    prev_node = g.predecessors(prev_node)[0]
                    if node_cnt >= flank:
                        break
                if node_cnt < flank: # and prev_node != st_node:
                    continue


            next_node = other_successor
            node_cnt = 0
            while g.in_degree(next_node) == 1 and g.out_degree(next_node) == 1:
                node_cnt += 1
                next_node = g.successors(next_node)[0]
                if node_cnt >= flank:
                    break

            if node_cnt < flank: # and next_node != first_node_of_repeat:
                continue

            rep = [first_node_of_repeat]
            next_node = first_node_of_repeat

            node_cnt = 0

            if g.in_degree(next_node) == 2 and g.out_degree(next_node) == 1:
                next_double_node = g.successors(next_node)[0]
                rep.append(next_double_node)
                prev_edge = g[next_node][next_double_node]
            else:
                next_double_node = next_node
                try:
                    assert not (g.in_degree(next_double_node) == 1 and g.out_degree(next_double_node) == 1)
                except:
                    print str(g.in_degree(next_node))
                    print str(g.out_degree(next_node))
                    raise

            while g.in_degree(next_double_node) == 1 and g.out_degree(next_double_node) == 1 and node_cnt < max_nodes:
                node_cnt += 1
                loop_len += abs(g[next_double_node][g.successors(next_double_node)[0]]['read_a_match_start'] - prev_edge['read_b_match_start'])
                next_double_node = g.successors(next_double_node)[0]
                rep.append(next_double_node)

            if print_debug:
                print "length in repeat " + str(loop_len-len_in_loop)

            if next_double_node == st_node and loop_len > MAX_PLASMID_LENGTH:
                if print_debug:
                    print 'success!'
                    print "length in loop " + str(loop_len)
                    print 'rep is:'
                    print rep
                    print 'in_node and other_successor:'
                    print in_node, other_successor
                resolve_rep(g,rep,in_node,other_successor)
    #             print next_double_node

                if node_cnt < 5:

                    tandem.append(rep)



                continue

    if len(tandem) > 0:
        with open('tandem.txt', 'w') as tandemout:
            for rep in tandem:
                tandemout.write(str(rep))


    return g




def y_pruning(G,flank):

    H = G.copy()

    y_nodes = set([x for x in H.nodes() if H.out_degree(x) > 1 and H.in_degree(x) == 1])

    pruned_count = 0

    for st_node in y_nodes:

        pruned = 0

        try:  
            H.predecessors(st_node)
        except:
            continue

        prev_node = H.predecessors(st_node)[0]

        node_cnt = 0
    
        while H.in_degree(prev_node) == 1 and H.out_degree(prev_node) == 1:
            node_cnt += 1
            prev_node = H.predecessors(prev_node)[0]
            if node_cnt >= flank:
                break
        if node_cnt < flank: # and prev_node != st_node:
            continue

        # if we got here, we probably have a Y, and not a collapsed repeat
        for vert in H.successors(st_node):
            if H.node[vert]['CFLAG'] == True:
                
                try:
                    H.remove_edge(st_node,vert)
                    H.remove_edge(rev_node(vert),rev_node(st_node))
                    pruned = 1
    
                except:
                    pass

        if pruned == 1:
            pruned_count += 1

    # print "Number of pruned Y's: "+str(pruned_count)


    return H


# In[72]:


def add_groundtruth(g,json_file,in_hinges,out_hinges):

    mapping = ujson.load(json_file)

    print 'getting mapping'
    mapped_nodes=0
    print str(len(mapping))
    print str(len(g.nodes()))

    slack = 500
    max_chr = 0

    chr_length_dict = {}

    for node in g.nodes():
        # print node
        node_base=node.split("_")[0]
        # print node_base

        #print node
        g.node[node]['normpos'] = 0
        if mapping.has_key(node_base):
            g.node[node]['chr'] = mapping[node_base][0][2]+1
            g.node[node]['aln_start'] = min (mapping[node_base][0][0],mapping[node_base][0][1])
            g.node[node]['aln_end'] = max(mapping[node_base][0][1],mapping[node_base][0][0])


            # max_chr = max(g.node[node]['chr'],max_chr)
            # mapped_nodes+=1
        else:
            # pass
            g.node[node]['chr'] =  0
            g.node[node]['aln_start'] = 1
            g.node[node]['aln_end'] = 1
#             g.node[node]['aln_strand'] = 0

        if node in in_hinges or node in out_hinges:
            g.node[node]['hinge'] = 1
        else:
            g.node[node]['hinge'] = 0

        if g.node[node]['chr'] in chr_length_dict:
            chr_length_dict[g.node[node]['chr']] = max(g.node[node]['aln_end'], chr_length_dict[g.node[node]['chr']])
        else:
            chr_length_dict[g.node[node]['chr']] = max(g.node[node]['aln_end'], 1)

    chr_list = sorted(chr_length_dict.items(), key=operator.itemgetter(1), reverse=True)

    max_chr_len1 = max([g.node[x]['aln_end'] for x in  g.nodes()])
    max_chr_multiplier = 10**len(str(max_chr_len1))
    print [x for x in chr_list]
    chr_set =[x [0] for x in chr_list]
    print chr_set
    # red_bk = 102
    # green_bk = 102
    # blue_bk = 102
    colour_list = ['red', 'lawngreen', 'deepskyblue', 'deeppink', 'darkorange', 'purple', 'gold', 'mediumblue',   'saddlebrown', 'darkgreen']
    for colour in colour_list:
        print  matplotlib.colors.colorConverter.to_rgb(colour)
    for index, chrom in enumerate(chr_set):
        node_set = set([x for x in  g.nodes() if g.node[x]['chr'] == chrom])
        print chrom


        max_chr_len = max([g.node[x]['aln_end'] for x in  g.nodes() if g.node[x]['chr'] == chrom])
        # max_chr_multiplier = 10**len(str(max_chr_len))


        if index < 10:
            rgb_tuple = matplotlib.colors.colorConverter.to_rgb(colour_list[index])
            red = int(255*rgb_tuple[0])
            green = int(255*rgb_tuple[1])
            blue = int(255*rgb_tuple[2])
        else:
            red = random.randint(0,255)
            # green = random.randint(0,255)
            blue = random.randint(0,255)
            brightness = 200
            green  = max(0,min( 255,brightness - int((0.2126 *red +  0.0722 *blue)/0.7152 )))

        red_bk = max(red-100,0)
        blue_bk = max(blue-100,0)
        green_bk = max(green-100,0)

        print red,blue,green
        for node in node_set:
            g.node[node]['normpos'] = g.node[node]['chr'] * max_chr_multiplier + (g.node[node]['aln_end']/float(max_chr_len))*max_chr_multiplier
            lamda = (g.node[node]['aln_end']/max_chr_len)
            nd_red = (1-lamda)*red + lamda*red_bk
            nd_green = (1-lamda)*green + lamda*green_bk
            nd_blue = (1-lamda)*blue + lamda*blue_bk
            g.node[node]['color'] = rgb2hex(nd_red, nd_green, nd_blue)
            g.node[node]['color_r'] = nd_red
            g.node[node]['color_g'] = nd_green
            g.node[node]['color_b'] = nd_blue

    # max_chr_len = len(str(max_chr))

    # div_num = float(10**(max_chr_len))

    # for node in g.nodes():
    #     g.node[node]['normpos'] = (g.node[node]['chr'] + g.node[node]['aln_end']/float(chr_length_dict[g.node[node]['chr']]))/div_num

    for edge in g.edges_iter():
        in_node=edge[0]
        out_node=edge[1]

#         if ((g.node[in_node]['aln_start'] < g.node[out_node]['aln_start'] and
#             g.node[out_node]['aln_start'] < g.node[in_node]['aln_end']) or
#             (g.node[in_node]['aln_start'] < g.node[out_node]['aln_end'] and
#             g.node[out_node]['aln_end'] < g.node[in_node]['aln_end'])):
#             g.edge[in_node][out_node]['false_positive']=0
#         else:
#             g.edge[in_node][out_node]['false_positive']=1


        if ((g.node[in_node]['aln_start'] < g.node[out_node]['aln_start'] and
            g.node[out_node]['aln_start'] < g.node[in_node]['aln_end']) or
            (g.node[in_node]['aln_start'] < g.node[out_node]['aln_end'] and
            g.node[out_node]['aln_end'] < g.node[in_node]['aln_end'])):
            g.edge[in_node][out_node]['false_positive']=0
        else:
            g.edge[in_node][out_node]['false_positive']=1

    return g


def mark_skipped_edges(G,skipped_name):

    with open (skipped_name) as f:
        for lines in f:
            lines1=lines.split()

            if len(lines1) < 5:
                continue

            e1 = (lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4])

            if e1 in G.edges():
                G.edge[lines1[0] + "_" + lines1[3]][lines1[1] + "_" + lines1[4]]['skipped'] = 1
                G.edge[lines1[1] + "_" + str(1-int(lines1[4]))][lines1[0] + "_" + str(1-int(lines1[3]))]['skipped'] = 1





def add_annotation(g,in_hinges,out_hinges):

    for node in g.nodes():

        if node in in_hinges:
            g.node[node]['hinge'] = 1
        elif node in out_hinges:
            g.node[node]['hinge'] = -1
        else:
            g.node[node]['hinge'] = 0

    return g




def add_chimera_flags(g,prefix):

    cov_flags = prefix + '.cov.flag'
    slf_flags = None

    for node in g.nodes():
        g.node[node]['CFLAG'] = False
    if slf_flags != None:
        g.node[node]['SFLAG'] = False

    node_set = set(g.nodes())
    num_bad_cov_reads = 0
    if cov_flags != None:
        with open(cov_flags,'r') as f:
            for line in f:
                node_name = line.strip()
                try: 
                    assert not ((node_name+'_0' in node_set and node_name+'_1' not in node_set)
                        or (node_name+'_0' not in node_set and node_name+'_1'  in node_set))
                except:
                    print node_name + ' is not symmetrically present in the graph input.'
                    raise
                if node_name+'_0' in node_set:
                    g.node[node_name+'_0']['CFLAG'] = True
                    g.node[node_name+'_1']['CFLAG'] = True
                    num_bad_cov_reads += 1
    print str(num_bad_cov_reads) + ' bad coverage reads.'

    num_bad_slf_reads = 0
    if slf_flags != None:
        with open(slf_flags,'r') as f:
            for line in f:
                node_name = line.strip()
                try: 
                    assert not ((node_name+'_0' in node_set and node_name+'_1' not in node_set)
                        or (node_name+'_0' not in node_set and node_name+'_1'  in node_set))
                except:
                    print node_name + ' is not symmetrically present in the graph input.'
                    raise
                if node_name+'_0' in node_set:
                    g.node[node_name+'_0']['SFLAG'] = True
                    g.node[node_name+'_1']['SFLAG'] = True
                    num_bad_slf_reads += 1
    print str(num_bad_slf_reads) + ' bad self aligned reads.'            






def connect_strands(g):

    for node in g.nodes():
        revnode = rev_node(node)
        g.add_edge(node,revnode)
        g.add_edge(revnode,node)

    return g





def create_bidirected(g):

    h = nx.DiGraph()

    for u in g.nodes():

        for successor in g.successors(u):

            tail_id, tail_orientation = u.split('_')
            head_id, head_orientation = successor.split('_')

            h.add_edge(tail_id,head_id,tail_or = int(tail_orientation),head_or = int(head_orientation),
                read_a_match_start=g.edge[u][successor]['read_a_match_start'],
                read_a_match_end=g.edge[u][successor]['read_a_match_end'],
                read_b_match_start=g.edge[u][successor]['read_b_match_start'],
                read_b_match_end=g.edge[u][successor]['read_b_match_end'])


    st_nodes = [x for x in g if g.in_degree(x) != 1 or g.out_degree(x) > 1]

    for st_node in st_nodes:

        for sec_node in g.successors(st_node):

            cur_node = st_node
            cur_id = cur_node.split('_')[0]
            next_node = sec_node
            next_id = next_node.split('_')[0]

            if next_id in h.successors(cur_id) and cur_id in h.successors(next_id):
                h.remove_edge(next_id,cur_id)

            while g.in_degree(next_node) == 1 and g.out_degree(next_node) == 1:

                cur_node = next_node
                cur_id = cur_node.split('_')[0]
                next_node = g.successors(next_node)[0]
                next_id = next_node.split('_')[0]
                # else:
                #     print 'not in h'

                if next_id in h.successors(cur_id) and cur_id in h.successors(next_id):
                    h.remove_edge(next_id,cur_id)
                else:
                    break


    return h




def create_bidirected2(g):

    h = nx.DiGraph()

    for u in g.nodes():

        for successor in g.successors(u):

            tail_id, tail_orientation = u.split('_')
            head_id, head_orientation = successor.split('_')

            h.add_edge(tail_id,head_id)

            # h.add_edge(tail_id,head_id,tail_or = int(tail_orientation),head_or = int(head_orientation),
            #     read_a_match_start=g.edge[u][successor]['read_a_match_start'],
            #     read_a_match_end=g.edge[u][successor]['read_a_match_end'],
            #     read_b_match_start=g.edge[u][successor]['read_b_match_start'],
            #     read_b_match_end=g.edge[u][successor]['read_b_match_end'])


    st_nodes = [x for x in g if g.in_degree(x) != 1 or g.out_degree(x) > 1]

    for st_node in st_nodes:

        for sec_node in g.successors(st_node):

            cur_node = st_node
            cur_id = cur_node.split('_')[0]
            next_node = sec_node
            next_id = next_node.split('_')[0]

            if next_id in h.successors(cur_id) and cur_id in h.successors(next_id):
                h.remove_edge(next_id,cur_id)

            while g.in_degree(next_node) == 1 and g.out_degree(next_node) == 1:

                cur_node = next_node
                cur_id = cur_node.split('_')[0]
                next_node = g.successors(next_node)[0]
                next_id = next_node.split('_')[0]
                # else:
                #     print 'not in h'

                if next_id in h.successors(cur_id) and cur_id in h.successors(next_id):
                    h.remove_edge(next_id,cur_id)
                else:
                    break


    return g

def update_progress(t):
    sys.stdout.write('\r{}%'.format(min(int(t),100)))
    sys.stdout.flush()

    
def update_fraction(flname,a,b):
    sys.stdout.write('\r[clip] Reading {}: {}/{}'.format(flname,a,b))
    sys.stdout.flush()
    

def write_graphml(g,prefix,suffix,suffix1):
    h = g.copy()
    connect_strands(h)
    nx.write_graphml(h, prefix+suffix+'.'+'suffix1'+'.graphml')


    
debug_mode = False

flname = sys.argv[1]
# flname = '../pb_data/ecoli_shortened/ecoli4/ecolii2.edges.hinges'

prefix = flname.split('.')[0]

hingesname = sys.argv[2]
# hingesname = '../pb_data/ecoli_shortened/ecoli4/ecolii2.hinge.list'


suffix = sys.argv[3]
DEL_TELOMERE = False
AGGRESSIVE_PRUNING = False

if len(sys.argv) >= 5:
    ini_file_path = sys.argv[4]
    config = configparser.ConfigParser()
    config.read(ini_file_path)
    try:
        MAX_PLASMID_LENGTH = config.getint('layout', 'max_plasmid_length')
        # print 'MAX_PLASMID_LENGTH in config '+str(MAX_PLASMID_LENGTH)
    except:
        MAX_PLASMID_LENGTH = 500000
        # print 'MAX_PLASMID_LENGTH '+str(MAX_PLASMID_LENGTH)
    try: 
        DEL_TELOMERE = config.getbool('layout','del_telomere')
    except:
        DEL_TELOMERE = False
    try: 
        AGGRESSIVE_PRUNING = config.getbool('layout','aggressive_pruning')
    except:
        AGGRESSIVE_PRUNING = False

else:
    MAX_PLASMID_LENGTH = 500000



if len(sys.argv)>=6:
    json_file = open(sys.argv[5])
else:
    json_file = None
# path = '../pb_data/ecoli_shortened/ecoli4/'
# suffix = 'i2'





# In[116]:

G = nx.DiGraph()

Ginfo = {}

print ("[clip] Started hinge clip ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
      
# print("Reading {}".format(flname))

file_length = 0
with open (flname) as f:
    for lines in f:
        file_length += 1

# print "{} edges total".format(file_length)        

# file_pc = file_length/100

file_count = 0

with open (flname) as f:
    for lines in f:
        lines1=lines.split()
        
        file_count+=1
        update_fraction(flname,file_count,file_length)

        if len(lines1) < 5:
            continue
       
            
#         if file_count % file_pc == 0:
#             update_progress(file_count*100/file_length)
#             print file_count*100/file_length
            
            
        e1 = (lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4])
        # print lines1
        # e1_match1 = abs(int(lines1[6].lstrip('['))-int(lines1[7].rstrip(']')))
        # e1_match2 = abs(int(lines1[8].lstrip('['))-int(lines1[9].rstrip(']')))
        e1_match_len = int(lines1[2])
        ra_match_start = int(lines1[6].lstrip('['))
        ra_match_end = int(lines1[7].rstrip(']'))
        rb_match_start = int(lines1[8].lstrip('['))
        rb_match_end = int(lines1[9].rstrip(']'))

        ra_match_start_raw = int(lines1[-4].lstrip('['))
        ra_match_end_raw = int(lines1[-3].rstrip(']'))
        rb_match_start_raw = int(lines1[-2].lstrip('['))
        rb_match_end_raw = int(lines1[-1].rstrip(']'))


        G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],
            hinge_edge=int(lines1[5]),length=e1_match_len,z=0,
            read_a_match_start=ra_match_start,read_a_match_end=ra_match_end,
            read_b_match_start=rb_match_start,read_b_match_end=rb_match_end,
            read_a_match_start_raw=ra_match_start_raw,read_a_match_end_raw=ra_match_end_raw,
            read_b_match_start_raw=rb_match_start_raw,read_b_match_end_raw=rb_match_end_raw)
        G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),
            hinge_edge=int(lines1[5]),length=e1_match_len,z=0,
            read_a_match_start=rb_match_start,read_a_match_end=rb_match_end,
            read_b_match_start=ra_match_start,read_b_match_end=ra_match_end,
            read_a_match_start_raw=rb_match_start_raw,read_a_match_end_raw=rb_match_end_raw,
            read_b_match_start_raw=ra_match_start_raw,read_b_match_end_raw=ra_match_end_raw)



        towrite = lines1[0] + "_" + lines1[3] +' '+ lines1[1] + "_" + lines1[4] +' '+ lines1[2]+' '+str(int(lines1[11][:-1])-int(lines1[10][1:]))+' '+str(int(lines1[13][:-1])-int(lines1[12][1:]))
        Ginfo[(lines1[0] + "_" + lines1[3],lines1[1] + "_" + lines1[4])] = towrite

        towrite= lines1[1] + "_" + str(1-int(lines1[4])) +' '+ lines1[0] + "_" + str(1-int(lines1[3])) +' '+ lines1[2]+' '+str(int(lines1[13][:-1])-int(lines1[12][1:]))+' '+str(int(lines1[11][:-1])-int(lines1[10][1:]))
        Ginfo[(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])))] = towrite
        
        
print("\n[clip] Graph with "+str(len(G))+" nodes built ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
    

if debug_mode:    
    print("Writing G00")
    nx.write_graphml(G, prefix+suffix+'.'+'G00'+'.graphml')


vertices=set()

in_hinges = set()
out_hinges = set()


file_length = 0
file_count=0
with open (hingesname) as f:
    for lines in f:
        file_length += 1

with open (hingesname) as f:
    for lines in f:
        file_count+=1
        update_fraction(flname,file_count,file_length)
        
        lines1=lines.split()

        if lines1[2] == '1':
            in_hinges.add(lines1[0]+'_0')
            out_hinges.add(lines1[0]+'_1')
        elif lines1[2] == '-1':
            in_hinges.add(lines1[0]+'_1')
            out_hinges.add(lines1[0]+'_0')
            
    print ""


add_annotation(G,in_hinges,out_hinges)

if debug_mode==True and os.path.isfile(prefix + '.cov.flag'):
    add_chimera_flags(G,prefix)

    
if debug_mode==True and os.path.isfile(flname.split('.')[0] + '.edges.skipped'):
    mark_skipped_edges(G,flname.split('.')[0] + '.edges.skipped')


if json_file!= None:
    add_groundtruth(G,json_file,in_hinges,out_hinges)

print ("[clip] Pruning of graph started ("+datetime.datetime.now().strftime("%H:%M:%S")+")") 
    
G0 = G.copy()

# Actual pruning, clipping and z deletion occurs below

# print(str(datetime.datetime.now()))      
# print("Clipping dead-ends...")  
print ("[clip] Clipping dead-ends ("+datetime.datetime.now().strftime("%H:%M:%S")+")") 

G0 = dead_end_clipping_sym(G0,10)

print("[clip] Number of nodes remaining: {}".format(len(G0)))  

# print(str(datetime.datetime.now()))      
# print("Clipping Z edges...")

print ("[clip] Clipping spurious edges ("+datetime.datetime.now().strftime("%H:%M:%S")+")") 

# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)
G1,G0 = z_clipping_sym(G0,6,set(),set())
# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)
# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)
# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)

print("[clip] Number of nodes remaining: {}".format(len(G1)))     

# print(str(datetime.datetime.now()))      
# print("Bursting bubbles...")
  
print ("[clip] Bursting bubbles ("+datetime.datetime.now().strftime("%H:%M:%S")+")")     
    
if DEL_TELOMERE:
    G1 = bubble_bursting_sym(G1,20)
    G1 = dead_end_clipping_sym(G1,20)
else:
    G1 = bubble_bursting_sym(G1,10)
    G1 = dead_end_clipping_sym(G1,5)
    
print("[clip] Number of nodes remaining: {}".format(len(G1)))  

           
print("[clip] Writing G0 ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
nx.write_graphml(G0, prefix+suffix+'.'+'G0'+'.graphml')

print("[clip] Writing G1 ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
nx.write_graphml(G1, prefix+suffix+'.'+'G1'+'.graphml')


G2 = G1.copy()

        
print("[clip] Condensing G1 ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
Gs = random_condensation_sym(G1,1000)


print("[clip] Performing loop resolution ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
loop_resolution(G2,500,50)

print("[clip] Writing G2 ("+datetime.datetime.now().strftime("%H:%M:%S")+")") 
nx.write_graphml(G2, prefix+suffix+'.'+'G2'+'.graphml')

   
print("[clip] Condensing G2 ("+datetime.datetime.now().strftime("%H:%M:%S")+")")    
G2s = random_condensation_sym(G2,1000)


nx.write_graphml(Gs, prefix+suffix+'.'+'Gs'+'.graphml')

nx.write_graphml(G2s, prefix+suffix+'.'+'G2s'+'.graphml')

print("[clip] Overlaying strands ("+datetime.datetime.now().strftime("%H:%M:%S")+")") 
Gc = connect_strands(Gs)
      

nx.write_graphml(Gc, prefix+suffix+'.'+'Gc'+'.graphml')

G2c = connect_strands(G2s)

nx.write_graphml(G2c, prefix+suffix+'.'+'G2c'+'.graphml')




if AGGRESSIVE_PRUNING:

    G3 = y_pruning(G2,10)

    G3 = dead_end_clipping_sym(G3,10)

    G3s = random_condensation_sym(G3,1000)

    G3c = connect_strands(G3s)

    nx.write_graphml(G3, prefix+suffix+'.'+'G2'+'.graphml')

    nx.write_graphml(G3s, prefix+suffix+'.'+'G3s'+'.graphml')

    nx.write_graphml(G3c, prefix+suffix+'.'+'G3c'+'.graphml')



print("[clip] Done ("+datetime.datetime.now().strftime("%H:%M:%S")+")") 




# G2b = create_bidirected2(G2)

# nx.write_graphml(G2b, prefix+suffix+'.'+'G2b'+'.graphml')




# H=prune_graph(G1,in_hinges,out_hinges)
# H=dead_end_clipping(H,5)


# I=prune_graph(H,in_hinges,out_hinges,True)
# I=dead_end_clipping(I,5)


# Gs = random_condensation(G1,2000)
# nx.write_graphml(Gs, path+'G'+suffix+'.graphml')
# write_graph(Gs,path+'G'+suffix+'.txt')

# Hs = random_condensation(H,2500)
# nx.write_graphml(Hs, path+'H'+suffix+'.graphml')
# write_graph(Hs,path+'H'+suffix+'.txt')

# Is = random_condensation(I,2500)
# nx.write_graphml(Is, path+'I'+suffix+'.graphml')
# write_graph(Is,path+'I'+suffix+'.txt')



