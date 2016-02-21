#!/usr/bin/python

# In[1]:

import networkx as nx
import sys

# In[2]:

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
        G.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4])
        G.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])))
        
        
nx.write_graphml(G, './out.graphml')

# In[ ]:



