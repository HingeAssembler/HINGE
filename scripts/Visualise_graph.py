#!/usr/bin/python

# In[1]:

import networkx as nx
import sys

# In[2]:
if len(sys.argv) >2:
    print "wrong usage.\n python Visualise_graph.py graph_edge_file [list_of_hinges]"

vertices=set()
with open (sys.argv[1]) as f:

    for lines in f:
        lines1=lines.split()
        #print lines1
        if len(lines1) < 5:
            continue
        #vertices.add(lines1[0])
        #vertices.add(str(lines1[1]))
        #vertices.add(str(lines1[0])+"_" + lines1[3])
        #vertices.add(str(lines1[1])+"_" + lines1[4])


# In[3]:

len(vertices)


# In[4]:

G = nx.DiGraph()
for vertex in vertices:
    G.add_node(vertex)


# In[5]:

with open (sys.argv[1]) as f:
    for lines in f:
        lines1=lines.split()
        print lines1
        if len(lines1) < 5:
            continue
        #print lines1
        G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4], hinge_edge=int(lines1[5]))
        G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))

try:        
    in_hinges = set()
    out_hinges = set()
    with open (sys.argv[2]) as f:

	for lines in f:
	    lines1=lines.split()

	    if lines1[2] == '1':
		in_hinges.add(lines1[0]+'_0')
		out_hinges.add(lines1[0]+'_1')
	    elif lines1[2] == '-1':
		in_hinges.add(lines1[0]+'_1')
		out_hinges.add(lines1[0]+'_0')  
            
    for node in G.nodes():
	if node in in_hinges and node in out_hinges:
	    G.node[node]['hinge']=100
	elif node in in_hinges:
	    G.node[node]['hinge']=10
	elif node in out_hinges:
	    G.node[node]['hinge']=-10
	else:
	    G.node[node]['hinge']=0
except:
    pass 	    
        
nx.write_graphml(G, './out.graphml')

# In[ ]:



