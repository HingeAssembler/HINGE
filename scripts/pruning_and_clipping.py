
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



# In[48]:

def merge_path(g,in_node,node,out_node):
        
    g.add_edge(in_node,out_node,hinge_edge = -1,false_positive = 0)        
    g.remove_node(node)

#     for edge in g.in_edges(in_node):
#         g.add_edge(edge[0],node,hinge_edge = -1)
    
#     for edge in g.out_edges(out_node):
#         g.add_edge(node,edge[1], hinge_edge = -1)
    
#     g.remove_node(in_node)
#     g.remove_node(out_node)
    


# In[121]:

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
   


flname = sys.argv[1]
# flname = '../pb_data/ecoli_shortened/ecoli4/ecolii2.edges.hinges'

prefix = flname.split('.')[0]

hingesname = sys.argv[2]
# hingesname = '../pb_data/ecoli_shortened/ecoli4/ecolii2.hinge.list'


json_file = open(sys.argv[3])


suffix = sys.argv[4]

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

        G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],hinge_edge=int(lines1[5]))
        G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))
        
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




# json_file = open('../pb_data/ecoli_shortened/ecoli4/ecoli.mapping.1.json')

add_groundtruth(G,json_file,in_hinges,out_hinges)


# In[ ]:

G0 = G.copy()

# Actual pruning, clipping and z deletion occurs below


G1=dead_end_clipping(G0,10)

G1=z_clipping(G1,5,in_hinges,out_hinges)

nx.write_graphml(G0, prefix+suffix+'.'+'G0'+'.graphml')

nx.write_graphml(G1, prefix+suffix+'.'+'G1'+'.graphml')


Gs = random_condensation(G1,2500)

nx.write_graphml(Gs, prefix+suffix+'.'+'Gs'+'.graphml')


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



