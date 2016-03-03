#!/usr/bin/python

import networkx as nx
import sys
from collections import Counter



def merge_simple_path(g):
    for node in g.nodes():
        #print g.in_degree(node), g.out_degree(node)
        if g.in_degree(node) == 1 and g.out_degree(node) == 1:
            
            in_node = g.in_edges(node)[0][0]
            out_node = g.out_edges(node)[0][1]
            if g.out_degree(in_node) == 1 and g.in_degree(out_node) == 1:
                if in_node != node and out_node != node and in_node != out_node:
                    merge_path(g,in_node,node,out_node)
                    
                              
def merge_two_nodes(g):
    for node in g.nodes():
        if g.in_degree(node) == 1 and g.out_degree(node) == 0:
            in_node = g.in_edges(node)[0][0]
            if g.out_degree(in_node) == 1:
                if in_node != node:
                    node_id = g.graph['aval']
                    g.graph['aval'] += 1
                    g.add_node(str(node_id), 
                        count = g.node[in_node]['count'] + g.node[node]['count'],
                        read = g.node[in_node]['read'] + ':' + g.node[node]['read'],
                        #aln_chr = g.node[node]['aln_chr']
                        )
                    g.remove_node(in_node)
                    g.remove_node(node)
            

def merge_path(g,in_node,node,out_node):
    #ov1 = find_overlap(g.node[in_node]['bases'], g.node[node]['bases'])
    #ov2 = find_overlap(g.node[node]['bases'], g.node[out_node]['bases'])
    
    node_id = g.graph['aval']
    g.graph['aval'] += 1
    #length = g.node[node]['length'] + g.node[in_node]['length'] + g.node[out_node]['length'] - ov1 - ov2
    #cov = (g.node[in_node]['cov'] * g.node[in_node]['length'] + g.node[node]['cov'] * g.node[node]['length']  + \
    #g.node[out_node]['cov'] * g.node[out_node]['length'])/float(length)
    #bases = g.node[in_node]['bases'][:-ov1] + g.node[node]['bases'] + g.node[out_node]['bases'][ov2:]
    
    g.add_node(str(node_id), 
        count = g.node[in_node]['count'] + g.node[node]['count'] + g.node[out_node]['count'],
        read = g.node[in_node]['read'] + ':' +  g.node[node]['read'] + ':' +g.node[out_node]['read'],
        #aln_chr = g.node[node]['aln_chr']
    )
    #g.add_node(str(node_id)+'-', bases = reverse_comp_bases(bases), length = length, cov = cov)
    #print g.node[str(node_id)]['chr']

    for edge in g.in_edges(in_node):
        g.add_edge(edge[0],str(node_id),st_pc=g.edge[edge[0]][edge[1]]['st_pc'],end_pc=g.edge[edge[0]][edge[1]]['end_pc'])

    
    for edge in g.out_edges(out_node):
        g.add_edge(str(node_id),edge[1],st_pc=g.edge[edge[0]][edge[1]]['st_pc'],end_pc=g.edge[edge[0]][edge[1]]['end_pc'])
    
        
    g.remove_node(in_node)
    g.remove_node(node)
    g.remove_node(out_node)
    
def input1(flname):
    g = nx.DiGraph()
    with open (flname) as f:
        for lines in f:
            lines1=lines.split()
            #print lines1
            if len(lines1) < 5:
                continue
            #print lines1
            g.add_edge(lines1[0] + "_" + lines1[3], lines1[1] + "_" + lines1[4], hinge_edge=int(lines1[5]))
            g.add_edge(lines1[1] + "_" + str(1-int(lines1[4])), lines1[0] + "_" + str(1-int(lines1[3])),hinge_edge=int(lines1[5]))
    return g
            
def input2(flname):
    g = nx.DiGraph()
    with open (flname) as f:
        for lines in f:
            lines1=lines.split()
            #print lines1
            g.add_edge(lines1[0], lines1[1])   
    return g

def run(filename, gt_file, n_iter):
    
    
    f=open(filename)
    line1=f.readline()
    print line1
    f.close()
    if len(line1.split()) !=2:
	   g=input1(filename)
    else:
	   g=input2(filename)
    
    read_to_chr_map={} 
    pos_dict = {}
    mapping_dict = {} 

    chr_lengths = {}
    for chr in range(14):
        chr_lengths[chr] = 1000

    with open(gt_file,'r') as f:
        for num, line in enumerate(f.readlines()):
            m = map(int, line.strip().split())
            # mapping_dict[num] = [min(m), max(m), int(m[0]>m[1])]
            read_to_chr_map[m[0]]= str(m[1])
            mapping_dict[num] = m[1]
            pos_dict[num] = [min(m[2],m[3]),max(m[2],m[3])]    
            # pos_dict[num] = [m[2],m[3],int(m[2]>m[3])]
            chr_lengths[m[1]] = max(chr_lengths[m[1]],max(m[2],m[3]))


    print nx.info(g)
    
    print "Chromosome lenghts:"
    print chr_lengths

    margin = 10000

    del_count = 0


    #print nx.info(g)
    print "Num reads read : "+str(len(read_to_chr_map))

    for cur_edge in g.edges():
        node0=int(cur_edge[0].split('_')[0])
        node1=int(cur_edge[1].split('_')[0])
        # g.edge[cur_edge[0]][cur_edge[1]]['st_pc'] = "{0:.2f}".format(1.0*pos_dict[node0][1]/chr_lengths[mapping_dict[node0]])
        # g.edge[cur_edge[0]][cur_edge[1]]['end_pc'] = "{0:.2f}".format(1.0*pos_dict[node1][0]/chr_lengths[mapping_dict[node1]])
        
        # st_pc is the "start percentage"; i.e., the percent location of edge[0] on its original chromosome
        # end_pc is the "end percentage"; i.e., the percent location of edge[1] on its original chromosome

        g.edge[cur_edge[0]][cur_edge[1]]['st_pc'] = 1.0*pos_dict[node0][1]/chr_lengths[mapping_dict[node0]]
        g.edge[cur_edge[0]][cur_edge[1]]['end_pc'] = 1.0*pos_dict[node1][0]/chr_lengths[mapping_dict[node1]]

    
    for node in g.nodes():
        nodeid=int(node.split('_')[0])

        if pos_dict[nodeid][0] < margin:
            g.remove_node(node)
            del_count += 1
            continue

        if pos_dict[nodeid][1] > chr_lengths[mapping_dict[nodeid]] - margin:
            g.remove_node(node)
            del_count += 1
            continue

        g.node[node]['count'] = 1
        g.node[node]['read'] = node
        #print str(nodeid), node,g.node[node]['chr']

    print "Deleted nodes: "+str(del_count)
        
        
    degree_sequence=sorted(g.degree().values(),reverse=True)
    print Counter(degree_sequence)
    for i in range(n_iter):
        for node in g.nodes():
            if g.in_degree(node) == 0:
                g.remove_node(node)
    
        print nx.info(g)
        degree_sequence=sorted(nx.degree(g).values(),reverse=True)
        print Counter(degree_sequence)
    
    degree_sequence=sorted(nx.degree(g).values(),reverse=True)
    print Counter(degree_sequence)
    
    
    g.graph['aval'] = 1000000000
    
    for i in range(5):
        merge_simple_path(g)
        degree_sequence=sorted(nx.degree(g).values(),reverse=True)
        print Counter(degree_sequence)
    
    h=nx.DiGraph()
    h.add_nodes_from(g)
    h.add_edges_from(g.edges())

    for cur_edge in h.edges():
        h.edge[cur_edge[0]][cur_edge[1]]['st_pc'] = g.edge[cur_edge[0]][cur_edge[1]]['st_pc']
        h.edge[cur_edge[0]][cur_edge[1]]['end_pc'] = g.edge[cur_edge[0]][cur_edge[1]]['end_pc']

    # h = g.copy()

    for node in g.nodes():
        reads_in_node=[int(x.split('_')[0]) for x in g.node[node]['read'].split(':')]
        try:
            chr_in_node=map(lambda x: read_to_chr_map[x], reads_in_node)
        except:
            print reads_in_node,g.node[node]['read']
            return
        chr_in_node_set=set(chr_in_node)
        if len(chr_in_node_set) ==1:
            h.node[node]['chr']=chr_in_node[0]
        else:
            h.node[node]['chr']= ':'.join(chr_in_node)

        h.node[node]['count']=g.node[node]['count']
        try:
            h.node[node]['read']=g.node[node]['read']
        except:
            pass


    try:
        import ujson
        mapping = ujson.load(open(filename.split('.')[0]+'.mapping.json'))
        
        print 'get mapping'
        
        for node in h.nodes():
            #print node
            if mapping.has_key(node):
                h.node[node]['aln_start'] = mapping[node][0]
                h.node[node]['aln_end'] = mapping[node][1]
                h.node[node]['aln_strand'] = mapping[node][2]
            else:
                h.node[node]['aln_start'] = 0
                h.node[node]['aln_end'] = 0
                h.node[node]['aln_strand'] = 0
                
    except:
        pass        


    
    nx.write_graphml(h, filename.split('.')[0]+'_condensed_annotated.graphml')
    nx.write_graphml(g, filename.split('.')[0]+'_G_condensed_annotated.graphml')
    
    print nx.number_weakly_connected_components(h)
    print nx.number_strongly_connected_components(h)
    
#

filename = sys.argv[1]
gt_file=sys.argv[2]
run(filename, gt_file,5)
