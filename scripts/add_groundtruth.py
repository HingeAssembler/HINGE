#!/usr/bin/python

import networkx as nx
import sys

graphml_file = sys.argv[1]
groundtruth_file = sys.argv[2]
graphml_file_w_groundtruth = sys.argv[3]
try:
    chromosome_to_consider= int(sys.argv[4])
except:
    chromosome_to_consider=None

g = nx.read_graphml(graphml_file)

print nx.info(g)

mapping_dict = {}

with open(groundtruth_file,'r') as f:
    for num, line in enumerate(f.readlines()):
        m = map(int, line.strip().split())
        # mapping_dict[num] = [min(m), max(m), int(m[0]>m[1])]
        mapping_dict[num] = [m[2],m[3],m[1]]
        
#print mapping_dict

max_len=0
for num in mapping_dict.keys():
    max_len=max(max_len,len(str(m[3])))


pow_mov=10**(max_len+1)
for node in g.nodes():
    #print node
    try:
        nodeid = int(node.split('_')[0])
        #print nodeid
        rev = int(node.split('_')[1])
        if chromosome_to_consider != None:
            g.node[node]['chromosome'] = 0
            if mapping_dict[nodeid][2]==chromosome_to_consider:
                g.node[node]['chromosome'] = mapping_dict[nodeid][2]+1
        else:
            g.node[node]['chromosome'] = mapping_dict[nodeid][2]+1
        
        if rev == 0:
            g.node[node]['aln_end'] =  mapping_dict[nodeid][2]*pow_mov+ mapping_dict[nodeid][1]
            g.node[node]['aln_start'] = mapping_dict[nodeid][2]*pow_mov + mapping_dict[nodeid][0]
            # g.node[node]['aln_strand'] = mapping_dict[nodeid][2]
        else:
            g.node[node]['aln_end'] = mapping_dict[nodeid][2]*pow_mov + mapping_dict[nodeid][1]
            g.node[node]['aln_start'] = mapping_dict[nodeid][2]*pow_mov+ mapping_dict[nodeid][0]
            # g.node[node]['aln_strand'] = 1-mapping_dict[nodeid][2]
            
    except:
        pass
        
nx.write_graphml(g, graphml_file_w_groundtruth)


