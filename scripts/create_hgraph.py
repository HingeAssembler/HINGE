#!/usr/bin/python

import networkx as nx
import random
import sys
from collections import Counter
import json



# This script creates a graphml file from the hgraph file


def read_graph(filename,gt_file):

    with open(gt_file) as f:
        read_dict=json.load(f)
    g = nx.DiGraph()

    with open (filename) as f:
        for lines in f:
            lines1=lines.split()
            g.add_node(lines1[0] + "_" + lines1[2])
            g.add_node(lines1[1] + "_" + lines1[3])
            if lines1[0] in read_dict:
                g.node[lines1[0] + "_" + lines1[2]]['aln_start']=min(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
                g.node[lines1[0] + "_" + lines1[2]]['aln_end']=max(read_dict[lines1[0]][0][0],read_dict[lines1[0]][0][1])
            else:
                g.node[lines1[0] + "_" + lines1[2]]['aln_start']=0
                g.node[lines1[0] + "_" + lines1[2]]['aln_end']=0
            if lines1[1] in read_dict:
                g.node[lines1[1] + "_" + lines1[3]]['aln_start']=min(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
                g.node[lines1[1] + "_" + lines1[3]]['aln_end']=max(read_dict[lines1[1]][0][0],read_dict[lines1[1]][0][1])
            else:
                g.node[lines1[1] + "_" + lines1[3]]['aln_start']=0
                g.node[lines1[1] + "_" + lines1[3]]['aln_end']=0

            g.node[lines1[0] + "_" + lines1[2]]['active']=1
            g.node[lines1[1] + "_" + lines1[3]]['active']=int(lines1[4])
            g.add_edge(lines1[0] + "_" + lines1[2], lines1[1] + "_" + lines1[3])
    
    
    nx.write_graphml(g, filename.split('.')[0]+'_hgraph.graphml')
    
    print nx.number_weakly_connected_components(g)
    print nx.number_strongly_connected_components(g)
  

if __name__ == "__main__":   

    read_graph(sys.argv[1],sys.argv[2])





