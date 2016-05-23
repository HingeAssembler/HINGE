
# coding: utf-8

# In[115]:

import networkx as nx
import random
import sys

import ujson



      

# print G.number_of_edges(),G.number_of_nodes()


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
#     H=nx.DiGraph()
    H = G.copy()
    start_nodes = set([x for x in H.nodes() if H.in_degree(x) ==0])
    
    for st_node in start_nodes:

        if st_node not in H.nodes():
            continue

        cur_path = [st_node]
        
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
                
        if print_debug:
            print '----2'
            print cur_path

            
        if len(cur_path) <= threshold:
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

    # g.add_edge(in_node,out_node,hinge_edge = -1,false_positive = 0)        

    if g.edge[in_node][node]['intersection'] == 1 and g.edge[node][out_node]['intersection'] == 1:
        g.add_edge(in_node,out_node,hinge_edge = -1,intersection = 1,z=0)
    else:
        g.add_edge(in_node,out_node,hinge_edge = -1,intersection = 0,z=0)

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

                        try:
                            merge_path(g,in_node,node,out_node)
                            merge_path(g,rev_node(out_node),rev_node(node),rev_node(in_node))
                        except:
                            pass
    
    if iter_cnt >= max_iter:
        print "couldn't finish sparsification"+str(len(g.nodes()))
                        
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
        read_a_start=g.edge[in_node][rep_path[0]]['read_a_start'],
        read_a_end=g.edge[in_node][rep_path[0]]['read_a_end'],
        read_b_start=g.edge[in_node][rep_path[0]]['read_b_start'],
        read_b_end=g.edge[in_node][rep_path[0]]['read_b_end']) 
    g.remove_edge(in_node,rep_path[0])

    g.add_edge(prefix+rep_path[-1],out_node,
        read_a_start=g.edge[rep_path[-1]][out_node]['read_a_start'],
        read_a_end=g.edge[rep_path[-1]][out_node]['read_a_end'],
        read_b_start=g.edge[rep_path[-1]][out_node]['read_b_start'],
        read_b_end=g.edge[rep_path[-1]][out_node]['read_b_end']) 
    g.remove_edge(rep_path[-1],out_node)
    

    g.add_edge(rev_node(prefix + rep_path[0]),rev_node(in_node),
        read_a_start=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_a_start'],
        read_a_end=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_a_end'],
        read_b_start=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_b_start'],
        read_b_end=g.edge[rev_node(rep_path[0])][rev_node(in_node)]['read_b_end']) 
    g.remove_edge(rev_node(rep_path[0]),rev_node(in_node))
    g.add_edge(rev_node(out_node),rev_node(prefix+rep_path[-1]),
        read_a_start=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_a_start'],
        read_a_end=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_a_end'],
        read_b_start=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_b_start'],
        read_b_end=g.edge[rev_node(out_node)][rev_node(rep_path[-1])]['read_b_end']) 
    g.remove_edge(rev_node(out_node),rev_node(rep_path[-1]))


   
    
    for i in range(0,len(rep_path)-1):
        g.add_edge(prefix+rep_path[i],prefix+rep_path[i+1],
            read_a_start=g.edge[rep_path[i]][rep_path[i+1]]['read_a_start'],
            read_a_end=g.edge[rep_path[i]][rep_path[i+1]]['read_a_end'],
            read_b_start=g.edge[rep_path[i]][rep_path[i+1]]['read_b_start'],
            read_b_end=g.edge[rep_path[i]][rep_path[i+1]]['read_b_end'])
        g.add_edge(rev_node(prefix+rep_path[i+1]),rev_node(prefix+rep_path[i]),
            read_a_start=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_a_start'],
            read_a_end=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_a_end'],
            read_b_start=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_b_start'],
            read_b_end=g.edge[rev_node(rep_path[i+1])][rev_node(rep_path[i])]['read_b_end'])
            
        


def loop_resolution(g,max_nodes,flank,print_debug = False):
    
    starting_nodes =  [x for x in g.nodes() if g.out_degree(x) == 2]

    if print_debug:
        print '----'
        print starting_nodes

    for st_node in starting_nodes:


        if g.out_degree(st_node) != 2:
            continue

        if print_debug:
            print '----'
            print st_node


        for first_node in g.successors(st_node):


            if g.out_degree(st_node) != 2:
                continue

            if print_debug:
                print '----'
                print first_node   


            other_successor = [x for x in g.successors(st_node) if x != first_node][0]
            next_node = other_successor
            node_cnt = 0
            while g.in_degree(next_node) == 1 and g.out_degree(next_node) == 1:
                node_cnt += 1
                next_node = g.successors(next_node)[0]
                if node_cnt >= flank:
                    break

            if node_cnt < flank:
                continue

            next_node = first_node
            if print_debug:
                print 'going on loop'

            node_cnt = 0
            while g.in_degree(next_node) == 1 and g.out_degree(next_node) == 1 and node_cnt < max_nodes:
                node_cnt += 1
                in_node = next_node
                next_node = g.successors(next_node)[0]

            if g.in_degree(next_node) == 2:
                prev_node = [x for x in g.predecessors(next_node) if x != in_node][0]
                node_cnt = 0
                while g.in_degree(prev_node) == 1 and g.out_degree(prev_node) == 1:
                    node_cnt += 1
                    prev_node = g.predecessors(prev_node)[0]
                    if node_cnt >= flank:
                        break
                if node_cnt < flank:
                    continue


            rep = [next_node] 

            node_cnt = 0

            if g.in_degree(next_node) == 2 and g.out_degree(next_node) == 1: 
                next_double_node = g.successors(next_node)[0]
                rep.append(next_double_node)
            else:
                next_double_node = next_node

            while g.in_degree(next_double_node) == 1 and g.out_degree(next_double_node) == 1 and node_cnt < max_nodes:
                node_cnt += 1
                next_double_node = g.successors(next_double_node)[0]
                rep.append(next_double_node)


            if next_double_node == st_node:
                if print_debug:
                    print 'success!'
                    print 'rep is:'
                    print rep
                    print 'in_node and other_successor:'
                    print in_node, other_successor
                resolve_rep(g,rep,in_node,other_successor)
    #             print next_double_node

                continue
        

    return g




# In[72]:


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
                read_a_start=g.edge[u][successor]['read_a_start'],
                read_a_end=g.edge[u][successor]['read_a_end'],
                read_b_start=g.edge[u][successor]['read_b_start'],
                read_b_end=g.edge[u][successor]['read_b_end'])


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
            #     read_a_start=g.edge[u][successor]['read_a_start'],
            #     read_a_end=g.edge[u][successor]['read_a_end'],
            #     read_b_start=g.edge[u][successor]['read_b_start'],
            #     read_b_end=g.edge[u][successor]['read_b_end'])


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




def write_graphml(g,prefix,suffix,suffix1):
    h = g.copy()
    connect_strands(h)
    nx.write_graphml(h, prefix+suffix+'.'+'suffix1'+'.graphml')



flname = sys.argv[1]
# flname = '../pb_data/ecoli_shortened/ecoli4/ecolii2.edges.hinges'

prefix = flname.split('.')[0]

hingesname = sys.argv[2]
# hingesname = '../pb_data/ecoli_shortened/ecoli4/ecolii2.hinge.list'


suffix = sys.argv[3]

if len(sys.argv)==5:
    json_file = open(sys.argv[4])
else:
    json_file = None
# path = '../pb_data/ecoli_shortened/ecoli4/'
# suffix = 'i2'



# In[116]:

G = nx.DiGraph()

Ginfo = {}

with open (flname) as f:
    for lines in f:
        lines1=lines.split()

        if len(lines1) < 5:
            continue

        
        e1 = (lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4])
        # print lines1
        # e1_match1 = abs(int(lines1[6].lstrip('['))-int(lines1[7].rstrip(']')))
        # e1_match2 = abs(int(lines1[8].lstrip('['))-int(lines1[9].rstrip(']')))
        e1_match_len = int(lines1[2])
        ra_start = int(lines1[6].lstrip('['))
        ra_end = int(lines1[7].rstrip(']'))
        rb_start = int(lines1[8].lstrip('['))
        rb_end = int(lines1[9].rstrip(']'))
        if e1 in G.edges():
            G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],
                hinge_edge=int(lines1[5]),intersection=1,length=e1_match_len,z=0,
                read_a_start=ra_start,read_a_end=ra_end,
                read_b_start=rb_start,read_b_end=rb_end)
            G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),
                hinge_edge=int(lines1[5]),intersection=1,length=e1_match_len,z=0,
                read_a_start=rb_start,read_a_end=rb_end,
                read_b_start=ra_start,read_b_end=ra_end)     
        else:
            G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],
                hinge_edge=int(lines1[5]),intersection=0,length=e1_match_len,z=0,
                read_a_start=ra_start,read_a_end=ra_end,
                read_b_start=rb_start,read_b_end=rb_end)
            G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),
                hinge_edge=int(lines1[5]),intersection=0,length=e1_match_len,z=0,
                read_a_start=rb_start,read_a_end=rb_end,
                read_b_start=ra_start,read_b_end=ra_end)




        
        towrite = lines1[0] + "_" + lines1[3] +' '+ lines1[1] + "_" + lines1[4] +' '+ lines1[2]+' '+str(int(lines1[11][:-1])-int(lines1[10][1:]))+' '+str(int(lines1[13][:-1])-int(lines1[12][1:]))
        Ginfo[(lines1[0] + "_" + lines1[3],lines1[1] + "_" + lines1[4])] = towrite
        
        towrite= lines1[1] + "_" + str(1-int(lines1[4])) +' '+ lines1[0] + "_" + str(1-int(lines1[3])) +' '+ lines1[2]+' '+str(int(lines1[13][:-1])-int(lines1[12][1:]))+' '+str(int(lines1[11][:-1])-int(lines1[10][1:]))
        Ginfo[(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])))] = towrite




vertices=set()

in_hinges = set()
out_hinges = set()

with open (hingesname) as f:

    for lines in f:
        lines1=lines.split()
       
        if lines1[2] == '1':
            in_hinges.add(lines1[0]+'_0')
            out_hinges.add(lines1[0]+'_1')
        elif lines1[2] == '-1':
            in_hinges.add(lines1[0]+'_1')
            out_hinges.add(lines1[0]+'_0')




add_annotation(G,in_hinges,out_hinges)

# try:
mark_skipped_edges(G,flname.split('.')[0] + '.edges.skipped')
# except:
#     print "some error here"
#     pass



# json_file = open('../pb_data/ecoli_shortened/ecoli4/ecoli.mapping.1.json')


if json_file!= None:
    add_groundtruth(G,json_file,in_hinges,out_hinges)


# In[ ]:

G0 = G.copy()

# Actual pruning, clipping and z deletion occurs below


G0 = dead_end_clipping_sym(G0,10)

# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)
G1,G0 = z_clipping_sym(G0,6,set(),set())
# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)
# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)
# G1=z_clipping_sym(G1,5,in_hinges,out_hinges)


G1 = bubble_bursting_sym(G1,10)

G1 = dead_end_clipping_sym(G1,5)

nx.write_graphml(G0, prefix+suffix+'.'+'G0'+'.graphml')
nx.write_graphml(G1, prefix+suffix+'.'+'G1'+'.graphml')


G2 = G1.copy()

Gs = random_condensation_sym(G1,1000)


loop_resolution(G2,500,50)

G2s = random_condensation_sym(G2,1000)





nx.write_graphml(G2, prefix+suffix+'.'+'G2'+'.graphml')

nx.write_graphml(Gs, prefix+suffix+'.'+'Gs'+'.graphml')

nx.write_graphml(G2s, prefix+suffix+'.'+'G2s'+'.graphml')

Gc = connect_strands(Gs)

nx.write_graphml(Gc, prefix+suffix+'.'+'Gc'+'.graphml')

G2c = connect_strands(G2s)

nx.write_graphml(G2c, prefix+suffix+'.'+'G2c'+'.graphml')

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



