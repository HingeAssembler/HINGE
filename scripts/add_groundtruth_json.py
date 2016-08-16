import networkx as nx
import sys
import json

graphml_file = sys.argv[1]
groundtruth_file = sys.argv[2]
graphml_file_w_groundtruth = sys.argv[3]

g = nx.read_graphml(graphml_file)

print nx.info(g)

with open(groundtruth_file) as f:
    read_dict=json.load(f)

max_len=0
for read  in read_dict:
    for aln_info in read_dict[read]:
        try:
            max_len=max(max_len,len(str(aln_info[0])))
            max_len=max(max_len,len(str(aln_info[1])))
        except:
            print 
            raise

pow_mov=10**(max_len+1)

for node in g.nodes():
    #print node
    nodeid = node.split('_')[0]
    #print nodeid
    rev = int(node.split('_')[1])
    if rev==1:
        nodeid+="'"

    if nodeid in read_dict:
        g.node[node]['chr'] = read_dict[nodeid][0][2]
        g.node[node]['aln_end'] =  pow_mov*read_dict[nodeid][0][2]+max(read_dict[nodeid][0][0],read_dict[nodeid][0][1])
#         g.node[node]['aln_start'] = pow_mov*read_dict[nodeid][0][2]+min(read_dict[nodeid][0][0],read_dict[nodeid][0][1])
#         g.node[node]['repeat']=0
#         if len (read_dict[nodeid]) >1 :
#             g.node[node]['repeat']=1
#             chrom_maps=set([aln[3] for aln in read_dict[nodeid]])
#             if len (chrom_maps) >  1:
#                 g.node[node]['repeat']=10
    else:
        g.node[node]['chr'] = -1
        g.node[node]['aln_end'] =  -1
#         g.node[node]['aln_start'] = -1
#         g.node[node]['repeat']=-1

nx.write_graphml(g, graphml_file_w_groundtruth)
