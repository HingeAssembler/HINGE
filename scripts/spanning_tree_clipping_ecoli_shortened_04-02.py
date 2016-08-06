
# coding: utf-8

# In[27]:

import networkx as nx
import random
import sys
# flname='../AwesomeAssembler/data/test_repeat/test_repeat.edges.hinges2'
flname = '../pb_data/ecoli_shortened/ecoli3/ecoli.edges.hinges'
# flname = '../AwesomeAssembler/data/test_triple_repeat/test_triple_repeat.edges.hinges.hinges.test'


# In[3]:

vertices=set()
with open (flname) as f:

    for lines in f:
        lines1=lines.split()
        #print lines1
        if len(lines1) < 5:
            continue


# In[4]:

G = nx.DiGraph()
for vertex in vertices:
    G.add_node(vertex)


# In[5]:

with open (flname) as f:
    for lines in f:
        lines1=lines.split()
        #print lines1
        #break
        if len(lines1) < 5:
            continue
        #print lines1
        G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4],hinge_edge=int(lines1[5]))
        G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))


# In[5]:

print G.number_of_edges(),G.number_of_nodes()


# In[6]:

node_list=[t for t in G.nodes_iter()]
high_dg_nodes=[node for node in node_list if (G.in_degree(node)>1 or G.out_degree(node)>1)]
high_in_dg_nodes=[node for node in node_list if (G.in_degree(node)>1) ]
high_out_dg_nodes=[node for node in node_list if (G.out_degree(node)>1) ]


# In[97]:

def write_graph(G,flname):
    with open(flname,'w') as f:
        for edge in G.edges_iter():
            f.write(str(edge[0])+'\t'+str(edge[1])+'\n')


# In[7]:

hinge_lines=[]
with open(flname) as f:
    for lines in f:
        lines1=lines.split()
        #print lines1
        if len(lines1) < 6:
            continue
        if lines1[5]=='1':
            hinge_lines.append(lines)


# In[8]:

hinge_set=set([line.split()[0]+'_0' for line in hinge_lines])
hinge_set=hinge_set.union([line.split()[0]+'_1' for line in hinge_lines])
hinge_set=hinge_set.union(set([line.split()[1]+'_1' for line in hinge_lines]))
hinge_set=hinge_set.union(set([line.split()[1]+'_0' for line in hinge_lines]))
print len(hinge_set)


# In[ ]:

import numpy as np
hinge_read=[x.split('_')[0] for x in in_hinges]
np.savetxt('../AwesomeAssembler/data/ecoli_normal7/used_hinges.txt', hinge_read,fmt='%s',delimiter='\n')


# In[15]:

hingesname = '../pb_data/ecoli_shortened/ecoli3/hinge_list.txt'
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


# In[12]:

in_hinges=hinge_set.intersection(in_hinges)
out_hinges=hinge_set.intersection(out_hinges)


# In[15]:

# hinges=['4842_0','4842_1', '4206_0', '4206_1','4717_0','4717_1', '1476_0', '1476_1','4844_0','4844_1', '3999_0', '3999_1']
# hinges=['12_0','12_1', '5068_0', '5068_1']
# in_hinges = ['12_0', '5068_0']
# out_hinges = ['12_1', '5068_1']

in_hinges = []
out_hinges = []


# In[17]:

def prune_graph(graph,in_hinges,out_hinges,reverse=False):
    H=nx.DiGraph()
    if reverse:
        G=nx.reverse(graph,copy=True)
    else:
        G=graph
    start_nodes = [x for x in G.nodes() if G.in_degree(x) ==0]
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
            print 'bad exit'
            return H
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


# In[18]:

def dead_end_clipping(G,threshold):
#     H=nx.DiGraph()
    H = G.copy()
    start_nodes = set([x for x in H.nodes() if H.in_degree(x) ==0])
    
    for st_node in start_nodes:
        cur_path = [st_node]
        
        if len(H.successors(st_node)) > 0:
            cur_node = H.successors(st_node)[0]
            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:
                cur_path.append(cur_node)
                cur_node = H.successors(cur_node)[0]
            
        if len(cur_path) <= threshold:
            for vertex in cur_path:
                H.remove_node(vertex)
    
    end_nodes = set([x for x in H.nodes() if H.out_degree(x) ==0])
    
    for end_node in end_nodes:
        cur_path = [end_node]
        if len(H.predecessors(end_node)) > 0:
            cur_node = H.predecessors(end_node)[0]
            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:
                cur_path.append(cur_node)
                cur_node = H.predecessors(cur_node)[0]
            
        if len(cur_path) <= threshold:
            for vertex in cur_path:
                H.remove_node(vertex)

    return H


# In[19]:

def z_clipping(G,threshold):
    H = G.copy()
    
    start_nodes = set([x for x in H.nodes() if H.out_degree(x) > 1])
    
    for st_node in start_nodes:
        for cur_node in H.successors(st_node):
            
            cur_path = [[st_node,cur_node]]

            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:
                cur_path.append([cur_node,H.successors(cur_node)[0]])
                cur_node = H.successors(cur_node)[0]
            
        if len(cur_path) <= threshold and H.in_degree(cur_node) > 1:
            for edge in cur_path:
                H.remove_edge(edge[0],edge[1])
            for j in range(len(cur_path)-1):
                H.remove_node(cur_path[j][1])
    
    end_nodes = set([x for x in H.nodes() if H.in_degree(x) > 1])
    
    for end_node in end_nodes:
        for cur_node in H.predecessors(end_node):
            
            cur_path = [[cur_node,end_node]]
        
            while H.in_degree(cur_node) == 1 and H.out_degree(cur_node) == 1:
                cur_path.append([H.predecessors(cur_node)[0],cur_node])
                cur_node = H.predecessors(cur_node)[0]
            
        if len(cur_path) <= threshold and H.out_degree(cur_node) > 1:
            for edge in cur_path:
                H.remove_edge(edge[0],edge[1])
            for j in range(len(cur_path)-1):
                H.remove_node(cur_path[j][0])

    return H


# In[41]:

def merge_path(g,in_node,node,out_node):
    #ov1 = find_overlap(g.node[in_node]['bases'], g.node[node]['bases'])
    #ov2 = find_overlap(g.node[node]['bases'], g.node[out_node]['bases'])
    
#     node_id = node
#     node_id = g.graph['aval']
#     g.graph['aval'] += 1
    #length = g.node[node]['length'] + g.node[in_node]['length'] + g.node[out_node]['length'] - ov1 - ov2
    #cov = (g.node[in_node]['cov'] * g.node[in_node]['length'] + g.node[node]['cov'] * g.node[node]['length']  + \
    #g.node[out_node]['cov'] * g.node[out_node]['length'])/float(length)
    #bases = g.node[in_node]['bases'][:-ov1] + g.node[node]['bases'] + g.node[out_node]['bases'][ov2:]
    
#     g.add_node(str(node_id))
    #g.add_node(str(node_id)+'-', bases = reverse_comp_bases(bases), length = length, cov = cov)
    
    for edge in g.in_edges(in_node):
        g.add_edge(edge[0],node,hinge_edge = -1)
    
    for edge in g.out_edges(out_node):
        g.add_edge(node,edge[1], hinge_edge = -1)
    
        
    g.remove_node(in_node)
#     g.remove_node(node)
    g.remove_node(out_node)
    


# In[42]:

def random_condensation(G,n_nodes):

    g = G.copy()
    
    while len(g.nodes()) > n_nodes:

        node = g.nodes()[random.randrange(len(g.nodes()))]

        if g.in_degree(node) == 1 and g.out_degree(node) == 1:

            # edge_1 = g.out_edges(node)[0]
            # edge_2 = g.in_edges(node)[0]

#             edge1 = g.out_edges(node)[0]
#             edge2 = g.in_edges(node)[0]

            # print g.edge[edge1[0]][edge1[1]]['hinge_edge']

#             if (g.edge[edge1[0]][edge1[1]]['hinge_edge'] == -1 and g.edge[edge2[0]][edge2[1]]['hinge_edge'] == -1):
            
            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    #print in_node, node, out_node
                    merge_path(g,in_node,node,out_node)
                        
                        
    return g


# In[94]:

#in_hinges_set = set(in_hinges)


# In[162]:

in_hinges_lst


# In[21]:

G=dead_end_clipping(G,10)
G=z_clipping(G,10)


# In[167]:

nx.write_graphml(G, './G8bout.graphml')
write_graph(G,'./G8out.txt')


# In[ ]:




# In[66]:

# in_hinges = ['38816_0']
# out_hinges = ['38816_1']
in_hinges = ['51084_0'] # arbitrary node
out_hinges = ['51084_1']


# In[22]:

in_hinges_lst = list(in_hinges.intersection(set(G.nodes())))
# set(G.nodes())
#out_hinges_set = set(out_hinges)
out_hinges_lst = list(out_hinges.intersection(set(G.nodes())))
H=prune_graph(G,in_hinges_lst,out_hinges_lst)
I=dead_end_clipping(H,5)


# In[40]:

nx.write_graphml(H, './H9out.graphml')
nx.write_graphml(I, './I9out.graphml')


# In[43]:


J = random_condensation(I,3000)


# In[44]:

nx.write_graphml(J, './Jout.graphml')


# In[103]:

len(out_hinges)


# In[107]:

J=prune_graph(I,in_hinges_lst,out_hinges_lst,True)
K=dead_end_clipping(J,10)


# In[108]:

nx.write_graphml(H, './H8out.graphml')
write_graph(H,'./H8out.txt')
nx.write_graphml(I, './I8out.graphml')
write_graph(I,'./I8out.txt')
# nx.write_graphml(J, './J7out.graphml')
# nx.write_graphml(K, './K7out.graphml')


# Debugging coverage

# In[121]:

hinge_read=set([int(x.split('_')[0]) for x in in_hinges])


# In[127]:

hinge_points={}
with open(hingesname) as f:
    for line in f:
        line1=line.split()
        if int(line1[0]) in hinge_read:
            hinge_points.setdefault(int(line1[0]),set()).add(int(line1[1]))


# In[128]:

hinge_points


# In[152]:

coverage_str=[]
cgs_str=[]
with open('/data/pacbio_assembly/AwesomeAssembler/data/ecoli_normal7/coverage.debug.txt', 'r') as f:
    for line_num, line in enumerate(f):
        if line_num in hinge_read:
            coverage_str.append(line)
with open('/data/pacbio_assembly/AwesomeAssembler/data/ecoli_normal7/coverage_gradient.debug.txt', 'r') as f:
    for line_num, line in enumerate(f):
        if line_num in hinge_read:
            cgs_str.append(line)


# In[153]:

get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt


for index, line in enumerate(coverage_str):
    read_loc=[]
    coverage=[]
    read_loc2=[]
    cgs=[]
    line1=line.strip().split()
    line2=cgs_str[index].strip().split()
    read_num=int(line1[0])
    max_cov=0
    for tup in line1[1:]:
        read_loc.append(int(tup.split(':')[0]))
        coverage.append(int(tup.split(':')[1]))
        max_cov=max(max_cov, int(tup.split(':')[1]))
    plt.subplot(211)
    plt.title(str(read_num+1))
    plt.plot(read_loc,coverage)
    plt.grid()
    for hinge in hinge_points[read_num]:
        plt.plot([hinge, hinge],[0, max_cov],linewidth =4)
    
    plt.subplot(212)
    for tup in line2[1:]:
        read_loc2.append(int(tup.split(':')[0]))
        cgs.append(int(tup.split(':')[1]))
    plt.plot(read_loc2,cgs)
    plt.grid()
    plt.show()


# In[ ]:



