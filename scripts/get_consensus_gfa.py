#!/usr/bin/env python

import sys
import os
import subprocess
from parse_read import *
import numpy as np
import networkx as nx
import itertools



filedir = sys.argv[1]
filename = sys.argv[2]
consensus_name = sys.argv[3]
in_graphml_name = filedir + '/' + filename +'_draft.graphml'


h = nx.read_graphml(in_graphml_name)


gfaname = filedir + '/' + filename +'_consensus.gfa'


consensus_contigs = []
try:
    with open(consensus_name) as f:
        for line in f:
            if line[0] != '>':
                consensus_contigs.append(line.strip())
except:
    pass

# for  i, vert in enumerate(h.nodes()):
#     print i, vert
#     try:
#         print i,len(h.node[vert]['path']), len(h.node[vert]['segment']), len(consensus_contigs[i])
#     except:
#         print len(h.nodes()), len(consensus_contigs)
#         raise

print 'Number of contigs'
print len(consensus_contigs), len(h.nodes())
with open(gfaname,'w') as f:
    f.write("H\tVN:Z:1.0\n")
    for j,vert in enumerate(h.nodes()):
        
        i = h.node[vert]['contig_id']
        print j, i

        seg = consensus_contigs[i]
        seg_line = "S\t"+vert+"\t"+seg + '\n'
        f.write(seg_line)
    for edge in h.edges():
        edge_line = "L\t"+edge[0]+"\t+\t"+edge[1]+"\t+\t0M\n"
        f.write(edge_line)

#last =  h.nodes()[-1]
#print h.node[last]
#path_last = h.node[last]['path']



#for i in range(len(path_last)-1):
#    read_a = path_last[i]
#    read_b = path_last[i+1]
#    print read_a, read_b, in_graph.edge[read_a][read_b]

# for i,node in enumerate(h.nodes()):
#      h.node[node]['path'] = ';'.join(h.node[node]['path'])
# nx.write_graphml(h,out_graphml_name)



