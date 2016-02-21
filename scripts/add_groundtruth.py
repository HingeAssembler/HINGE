#!/usr/bin/python

import networkx as nx
import sys

graphml_file = sys.argv[1]
groundtruth_file = sys.argv[2]
graphml_file_w_groundtruth = sys.argv[3]

g = nx.read_graphml(graphml_file)

print nx.info(g)

mapping_dict = {}

with open(groundtruth_file,'r') as f:
    for num, line in enumerate(f.readlines()):
        m = map(int, line.strip().split())
        mapping_dict[num] = [min(m), max(m), int(m[0]>m[1])]
        
#print mapping_dict

for node in g.nodes():
    #print node
    try:
        nodeid = int(node.split('_')[0])
        #print nodeid
        rev = int(node.split('_')[1])
        
        if rev == 0:
            g.node[node]['aln_start'] = mapping_dict[nodeid][0]
            g.node[node]['aln_end'] = mapping_dict[nodeid][1]
            g.node[node]['aln_strand'] = mapping_dict[nodeid][2]
        else:
            g.node[node]['aln_start'] = mapping_dict[nodeid][0]
            g.node[node]['aln_end'] = mapping_dict[nodeid][1]
            g.node[node]['aln_strand'] = 1-mapping_dict[nodeid][2]
            
    except:
        pass
        
nx.write_graphml(g, graphml_file_w_groundtruth)