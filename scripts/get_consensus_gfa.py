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
map_filename = filedir + '/draft_map.txt'


g = nx.read_graphml(in_graphml_name)


gfaname = filedir + '/' + filename +'_consensus.gfa'
cols = np.loadtxt(map_filename, dtype=str,usecols=(1,))
del_contigs = np.nonzero(cols == 'Deleted')[0]


# consensus_contigs = []
# i = 0
# try:
#     with open(consensus_name) as f:
#         for  line in f:
#             if line[0] != '>':
#                 consensus_contigs.append(line.strip())
#                 i += 1
#                 while i in set(del_contigs):
#                     consensus_contigs.append('')
#                     print len()
#                     i += 1
# except:
#     pass

del_contig_ptr = 0
cols = np.loadtxt('/data/pacbio_assembly/pb_data/yeast/draft_map.txt', dtype=str,usecols=(1,))
del_contigs = np.nonzero(cols == 'Deleted')[0]

consensus_contigs = []
i = 0
with open(consensus_name) as f:
    for  line in f:
        if line[0] != '>':
            while del_contig_ptr < len(del_contigs) :
                if len(consensus_contigs) == del_contigs[del_contig_ptr]:
                    consensus_contigs.append('')
                    del_contig_ptr += 1
                else:
                    break
            consensus_contigs.append(line.strip())
            i += 1


nodes_to_keep = [x for x in g.nodes() if consensus_contigs[g.node[x]['contig_id']] != '' ]
h = g.subgraph(nodes_to_keep)


# for  i, vert in enumerate(h.nodes()):
#     print i, vert
#     try:
#         print i,len(h.node[vert]['path']), len(h.node[vert]['segment']), len(consensus_contigs[i])
#     except:
#         print len(h.nodes()), len(consensus_contigs)
#         raise

print 'Number of contigs'
print len(consensus_contigs), len(h.nodes())
# print [len(x) for x in consensus_contigs]


with open(gfaname,'w') as f:
    f.write("H\tVN:Z:1.0\n")
    for j,vert in enumerate(h.nodes()):
        
        i = h.node[vert]['contig_id']
        # print j, i

        seg = consensus_contigs[i]
        print(len(seg))
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



