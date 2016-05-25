#!/usr/bin/python
import networkx as nx
import sys

filename = sys.argv[1]
outfile = filename.split('.')[0] + ".edges.list"

g = nx.read_graphml(filename)
print nx.info(g)

def get_unitig(g, node):
    start_node = node
    path = [start_node]
    path_before = []
    path_after = []
    
    node = start_node
    while g.in_degree(node) == 1:
        node = g.in_edges(node)[0][0]
        if g.out_degree(node) == 1:
            if node in path_before:
                break
            else:
                path_before.append(node)
        else:
            break
            
    node = start_node
    while g.out_degree(node) == 1:
        node = g.out_edges(node)[0][1]
        if g.in_degree(node) == 1:
            if node in path_after or node in path_before:
                break
            else:
                path_after.append(node)
        else:
            break
        
    path = path_before[::-1] + path + path_after
    
    return path
    
def get_unitigs(g):
    paths = []
    node_set = set(g.nodes())
    while len(node_set) > 0:
        node = next(iter(node_set))
        path = get_unitig(g, node)
        node_set -= set(path)
        paths.append(path)
    return paths
        
        
paths = get_unitigs(g)

print len(paths)

with open(outfile, 'w') as f:
    for i,path in enumerate(paths):
        f.write('>Unitig%d\n'%(i))
        for j in range(len(path)-1):
            d =  g.get_edge_data(path[j],path[j+1])
            f.write('%s %s %s %s %d %d %d %d %d\n'%((path[j].split('_'))[0],path[j].split('_')[1]  , path[j+1].split('_')[0], path[j+1].split('_')[1], -d['read_a_start_raw'] + d['read_a_end_raw'] - d['read_b_start_raw'] + d['read_b_end_raw'], d['read_a_start_raw'], d['read_a_end_raw'], d['read_b_start_raw'], d['read_b_end_raw']))
    f.close()