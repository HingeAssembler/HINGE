#!/usr/bin/python

import networkx as nx
import random
import sys
from collections import Counter
import ujson



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



# In[9]:

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





def merge_path(g,in_node,node,out_node):
        
    g.add_edge(in_node,out_node,hinge_edge = -1,false_positive = 0)        
    g.remove_node(node)



def merge_a_to_b(g,node_a,node_b):
        
    if node_a not in g.nodes() or node_b not in g.nodes():
        return

    for node in g.predecessors(node_a):
        if node != node_b:
            g.add_edge(node,node_b,hinge_edge = 1,false_positive = 0)

    for node in g.successors(node_a):
        if node != node_b:
            g.add_edge(node_b,node,hinge_edge = 1,false_positive = 0)
   
    g.remove_node(node_a)


def random_condensation(G,n_nodes):

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



def add_groundtruth(g,json_file,in_hinges,out_hinges):
    
    mapping = ujson.load(json_file)

    print 'getting mapping'
    mapped_nodes=0
    print str(len(mapping)) 
    print str(len(g.nodes()))
    
    slack = 500
           
            
    for node in g.nodes():
        # print node
        node_base=node.split("_")[0]
        # print node_base

        #print node
        if mapping.has_key(node_base):
            g.node[node]['aln_start'] = min (mapping[node_base][0][0],mapping[node_base][0][1])
            g.node[node]['aln_end'] = max(mapping[node_base][0][1],mapping[node_base][0][0])
#             g.node[node]['chr'] = mapping[node_base][0][2]
            mapped_nodes+=1
        else:
            # pass
            g.node[node]['aln_start'] = 0
            g.node[node]['aln_end'] = 0
#             g.node[node]['aln_strand'] = 0
            
        if node in in_hinges or node in out_hinges:
            g.node[node]['hinge'] = 1
        else:
            g.node[node]['hinge'] = 0


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
   





def read_graph(edges_file,hg_file,gt_file, hinge_file):


    prefix = edges_file.split('.')[0]

    with open(gt_file) as f:
        read_dict = ujson.load(f)

    g = nx.DiGraph()

    hinge_nodes = []
    hinge_pos = {}


    with open (hinge_file) as f:
        for lines in f:
            lines1=lines.split()
            hinge_nodes.append(lines1[0] + "_0_" + lines1[1])
            hinge_nodes.append(lines1[0] + "_1_" + lines1[1])
            

            # if lines1[0] not in hinge_pos:
            #     hinge_pos[lines]
            # hinge_pos[lines1[0]] = lines1[1]


    # with open (hg_file) as f:
    #     for lines in f:
    #         lines1=lines.split()
    #         g.add_node(lines1[0] + "_" + lines1[2])
    #         g.add_node(lines1[1] + "_" + lines1[3])
    #         if lines1[0] in read_dict:
    #             g.node[lines1[0] + "_" + lines1[2]]['aln_start']=min(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
    #             g.node[lines1[0] + "_" + lines1[2]]['aln_end']=max(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
    #         else:
    #             g.node[lines1[0] + "_" + lines1[2]]['aln_start']=0
    #             g.node[lines1[0] + "_" + lines1[2]]['aln_end']=0
    #         if lines1[1] in read_dict:
    #             g.node[lines1[1] + "_" + lines1[3]]['aln_start']=min(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
    #             g.node[lines1[1] + "_" + lines1[3]]['aln_end']=max(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
    #         else:
    #             g.node[lines1[1] + "_" + lines1[3]]['aln_start']=0
    #             g.node[lines1[1] + "_" + lines1[3]]['aln_end']=0

    #         if lines1[0] in hinge_nodes:
    #             g.node[lines1[0] + "_" + lines1[2]]['active']=2
    #         else:
    #             g.node[lines1[0] + "_" + lines1[2]]['active']=1

    #         if lines1[1] in hinge_nodes:
    #             g.node[lines1[1] + "_" + lines1[3]]['active']=2
    #         else:
    #             g.node[lines1[1] + "_" + lines1[3]]['active']=int(lines1[4])


    #         g.add_edge(lines1[0] + "_" + lines1[2], lines1[1] + "_" + lines1[3], rev = int(lines1[5]))


    # need to construct double stranded hinge graph, so that proper mapping can be found

    with open (hg_file) as f:
        for lines in f:
            lines1=lines.split()

            nodeA0 = lines1[0] + "_0_"+ lines1[2]
            nodeA1 = lines1[0] + "_1_"+ lines1[2]
            nodeB0 = lines1[1] + "_0_"+ lines1[3]
            nodeB1 = lines1[1] + "_1_"+ lines1[3]            

            nodeA0short = lines1[0] + "_0"
            nodeA1short = lines1[0] + "_1"
            nodeB0short = lines1[1] + "_0"
            nodeB1short = lines1[1] + "_1"

            g.add_node(nodeA0)
            g.add_node(nodeA1)
            g.add_node(nodeB0)
            g.add_node(nodeB1)

            if nodeA0short not in hinge_pos:
                hinge_pos[nodeA0short] = [int(lines1[2])]
                hinge_pos[nodeA1short] = [int(lines1[2])]
            elif lines1[2] not in hinge_pos[nodeA0short]:
                hinge_pos[nodeA0short].append(int(lines1[2]))
                hinge_pos[nodeA1short].append(int(lines1[2]))

            if nodeB0 not in hinge_pos:
                hinge_pos[nodeB0short] = [int(lines1[3])]
                hinge_pos[nodeB1short] = [int(lines1[3])]
            elif lines1[3] not in hinge_pos[nodeB0short]:
                hinge_pos[nodeB0short].append(int(lines1[3]))
                hinge_pos[nodeB1short].append(int(lines1[3]))


            if lines1[0] in read_dict:
                g.node[lines1[0] + "_0_"+ lines1[2]]['aln_start']=min(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
                g.node[lines1[0] + "_0_"+ lines1[2]]['aln_end']=max(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
                g.node[lines1[0] + "_1_"+ lines1[2]]['aln_start']=min(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
                g.node[lines1[0] + "_1_"+ lines1[2]]['aln_end']=max(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
            else:
                g.node[lines1[0] + "_0_"+ lines1[2]]['aln_start']=0
                g.node[lines1[0] + "_0_"+ lines1[2]]['aln_end']=0
                g.node[lines1[0] + "_1_"+ lines1[2]]['aln_start']=0
                g.node[lines1[0] + "_1_"+ lines1[2]]['aln_end']=0
            if lines1[1] in read_dict:
                g.node[lines1[1] + "_0_"+ lines1[3]]['aln_start']=min(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
                g.node[lines1[1] + "_0_"+ lines1[3]]['aln_end']=max(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
                g.node[lines1[1] + "_1_"+ lines1[3]]['aln_start']=min(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
                g.node[lines1[1] + "_1_"+ lines1[3]]['aln_end']=max(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
            else:
                g.node[lines1[1] + "_0_"+ lines1[3]]['aln_start']=0
                g.node[lines1[1] + "_0_"+ lines1[3]]['aln_end']=0
                g.node[lines1[1] + "_1_"+ lines1[3]]['aln_start']=0
                g.node[lines1[1] + "_1_"+ lines1[3]]['aln_end']=0

            if nodeA0 in hinge_nodes:
                g.node[nodeA0]['active']=2
                g.node[nodeA1]['active']=2
            else:
                g.node[nodeA0]['active']=1
                g.node[nodeA1]['active']=1

            if nodeB0 in hinge_nodes:
                g.node[nodeB0]['active']=2
                g.node[nodeB1]['active']=2
            else:
                g.node[nodeB0]['active']=int(lines1[4])
                g.node[nodeB1]['active']=int(lines1[4])


            if int(lines1[5]) == 1: # reverse match
                g.add_edge(nodeA0, nodeB1)
                g.add_edge(nodeA1, nodeB0)
            else:
                g.add_edge(nodeA0,nodeB0)
                g.add_edge(nodeA1,nodeB1)

    
    
    # nx.write_graphml(g, filename.split('.')[0]+'_hgraph.graphml')
    
    # for c in nx.connected_components(g):
        # print len(c)

    hinge_mapping = {}

    for c in nx.weakly_connected_components(g):

        if len(c) > 10:

            component_sink = -1

            for node in c:
                if g.out_degree(node) == 0 and g.node[node]['active']== 2 and component_sink == -1:
                    component_sink = node
                elif g.out_degree(node) == 0 and g.node[node]['active']== 2 and g.in_degree(node) > g.in_degree(component_sink):
                    component_sink = node

            if component_sink != -1:
                g.node[component_sink]['active'] = 3
            else:
                component_sink = list(c)[0]

                # sink_shortname = component_sink.split('_')[0]

            for node in c:
                hinge_mapping[node] = component_sink

        else:

            for node in c:
                g.node[node]['active'] = -1


    nx.write_graphml(g, hg_file.split('.')[0]+'_hgraph2.graphml')


    # print nx.number_weakly_connected_components(g)
    # print nx.number_strongly_connected_components(g)


    G = nx.DiGraph()

    merging = 1

    if merging == 0:

        with open (edges_file) as f:
            for lines in f:
                lines1=lines.split()

                if len(lines1) < 6:
                    continue

                if int(lines1[5]) != 0:

                    if int(lines1[5]) == 1: 
                        # nodeB_id = lines1[1]+"_"+lines1[4]
                        # hingepos = int(lines1[6])
                        nodeA_id = lines1[0] + "_" + lines1[3]

                        hinge_node = lines1[1]+"_"+lines1[4] + '_' + lines1[6]

                        print hinge_node

                        eff_hinge = hinge_mapping[hinge_node]

                        eff_b = eff_hinge.split('_')

                        if eff_b[0] + "_" + eff_b[1] != lines1[0] + "_" + lines1[3]:
                            G.add_edge(lines1[0] + "_" + lines1[3], eff_b[0] + "_" + eff_b[1],hinge_edge=1)
                            G.add_edge(eff_b[0] + "_" + str(1-int(eff_b[1])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=1)
                        else:

                            G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],hinge_edge=1)
                            G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=1)


                    elif int(lines1[5]) == -1:

                        hinge_node = lines1[0]+"_"+lines1[3] + '_' + lines1[6]

                        eff_hinge = hinge_mapping[hinge_node]

                        eff_b = eff_hinge.split('_')

                        if eff_b[0] + "_" + eff_b[1] != lines1[1] + "_" + lines1[4] :
                            G.add_edge(eff_b[0] + "_" + eff_b[1], lines1[1] + "_" + lines1[4] ,hinge_edge=1)
                            G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), eff_b[0] + "_" + str(1-int(eff_b[1])),hinge_edge=1)
                        else:
                            G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],hinge_edge=1)
                            G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=1)


                    # if nodeA_id in hinge_pos:
                    #     print "found A"
                    # else:
                    #     print "didnt find A"


                    # if nodeB_id in hinge_pos:

                    #     print "Node B of hinged match IS in hinge_pos"

                    #     hinge_found = False

                    #     hingepos = hinge_pos[nodeB_id][0]
                    #     for candidate_pos in hinge_pos[nodeB_id]:
                    #         if (abs(candidate_pos - int(lines1[8][1:])) < 200):
                    #             hingepos = candidate_pos

                    #             print "Matching hinge found"

                    #             hinge_found = True

                    #     if not hinge_found:
                    #         print "not found"
                    #         print lines1[8][1:], hinge_pos[nodeB_id]


                    #     hinge_node = nodeB_id + '_' + str(hingepos)
                    #     eff_hinge = hinge_mapping[hinge_node]

                    #     eff_b = eff_hinge.split('_')

                    #     G.add_edge(lines1[0] + "_" + lines1[3], eff_b[0] + "_" + eff_b[1],hinge_edge=1)
                    #     G.add_edge(eff_b[0] + "_" + str(1-int(eff_b[1])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=1)


                else:

                    G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],hinge_edge=int(lines1[5]))
                    G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))



    else:

        to_be_merged = []

        with open (edges_file) as f:
            for lines in f:
                lines1=lines.split()

                if len(lines1) < 6:
                    continue

                G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],hinge_edge=int(lines1[5]))
                G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))


                if int(lines1[5]) != 0:

                    if int(lines1[5]) == 1: 

                        to_be_merged.append([lines1[1],lines1[6]])

                    elif int(lines1[5]) == -1:

                        to_be_merged.append([lines1[0],lines1[6]])

        for pair in to_be_merged:

            sink_node_long = hinge_mapping[pair[0]+'_0_'+pair[1]]
            sink_node = sink_node_long.split('_')[0]+'_'+sink_node_long.split('_')[1]

            if pair[0]+'_0' != sink_node:
                merge_a_to_b(G,pair[0]+'_0',sink_node)

            sink_node_long = hinge_mapping[pair[0]+'_1_'+pair[1]]
            sink_node = sink_node_long.split('_')[0]+'_'+sink_node_long.split('_')[1]

            if pair[0]+'_1' != sink_node:
                merge_a_to_b(G,pair[0]+'_1',sink_node)
            

    in_hinges = set()
    out_hinges = set()

    with open (hinge_file) as f:

        for lines in f:
            lines1=lines.split()
           
            if lines1[2] == '1':
                in_hinges.add(lines1[0]+'_0')
                out_hinges.add(lines1[0]+'_1')
            elif lines1[2] == '-1':
                in_hinges.add(lines1[0]+'_1')
                out_hinges.add(lines1[0]+'_0')


    json_file = open(gt_file)
    add_groundtruth(G,json_file,in_hinges,out_hinges)


    G0 = G.copy()

    nx.write_graphml(G0, prefix+'.'+'G0_merged'+'.graphml')


    G0s = random_condensation(G0,3500)

    nx.write_graphml(G0s, prefix+'.'+'G0s_merged'+'.graphml')


    G1=dead_end_clipping(G0,10)

    G1=z_clipping(G1,5,in_hinges,out_hinges)

    nx.write_graphml(G1, prefix+'.'+'G1_merged'+'.graphml')


    Gs = random_condensation(G1,2500)

    nx.write_graphml(Gs, prefix+'.'+'Gs_merged'+'.graphml')

  

if __name__ == "__main__":   

    read_graph(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])





