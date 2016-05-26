#!/usr/bin/python
import networkx as nx
import sys
import itertools

filename = sys.argv[1]
outfile = filename.split('.')[0] + ".edges.list"

g = nx.read_graphml(filename)
print nx.info(g)


def get_circle(g,node,vertices_of_interest):
    cur_path = [node]
    cur_vertex = g.successors(node)[0]

    i = 0
    while cur_vertex != node:
        cur_path.append(cur_vertex)
        try:
            assert len(g.successors(cur_vertex)) == 1
        except:
            print g.successors(cur_vertex), cur_vertex, node
            print cur_vertex in vertices_of_interest
            raise
        successor = g.successors(cur_vertex)[0]
        cur_vertex = successor
        
    cur_path.append(cur_vertex)
    


    return cur_path

    
def get_unitigs(g):
    paths = []
    num_paths = 0
    node_set = set(g.nodes())

    
    vertices_of_interest = set([x for x in g if g.in_degree(x) != 1 or g.out_degree(x) != 1])
    vertices_used = set(vertices_of_interest)
    for start_vertex in vertices_of_interest:
        first_out_vertices = g.successors(start_vertex)
        print first_out_vertices
        for vertex in first_out_vertices:
            cur_path = [start_vertex]
            cur_vertex = vertex
            while cur_vertex not in vertices_of_interest:
                successor = g.successors(cur_vertex)[0]
                cur_path.append(cur_vertex)
                predecessor = cur_vertex
                cur_vertex = successor
                
            cur_path.append(cur_vertex)
            vertices_used = vertices_used.union(set(cur_path))
            paths.append(cur_path)

    print len(node_set)
    print len(vertices_used)

    while len(node_set-vertices_used) > 0:
        node = list(node_set-vertices_used)[0]
        # print list(node_set-vertices_used)
        # # print vertices_of_interest
        # # print len(node_set-vertices_used)
        # break
        path = get_circle(g, node, vertices_of_interest)
        vertices_used = vertices_used.union(set(path))
        if len(path) > 1:
            paths.append(path)
    print len(paths) 
    # print paths
    print "paths"
    return paths
        
        
paths = get_unitigs(g)

print len(paths)

h = nx.DiGraph()
for i, path in enumerate(paths):
    h.add_node(i)
    h.node[i]['path'] = path

vertices_of_interest = set([x for x in g if g.in_degree(x) != 1 or g.out_degree(x) != 1])
for vertex in vertices_of_interest:
    successors = [x for x in h.nodes() if h.node[x]['path'][0] == vertex]
    predecessors = [x for x in h.nodes() if h.node[x]['path'][-1] == vertex]
    print successors,predecessors
    assert len(successors)==1 or len(predecessors)==1
    for succ, pred in itertools.product(successors,predecessors):
        h.add_edge(pred,succ)

    # if vertex.split('_') == '0':
    #     if len(predecessors) == 1:
    #         for succ in successors:
    #             rel_suc = h.node[succ]['path'][1]
    #             d = g.get_edge_data(vertex,rel_suc)
    #             h.edge[predecessors[0]][succ]['start_pos'] = d['read_a_end']
    #             h.edge[predecessors[0]][succ]['weight'] = d['read_a_end_raw'] - d['read_a_start']
    #     if len(successors) == 1:
    #         for pred in predecessors:
    #             rel_pred = h.node[pred]['path'][-2]
    #             d = g.get_edge_data(rel_pred,vertex)
    #             h.edge[predecessors[0]][succ][''] = d['read_a_end_raw'] - d['read_a_start']


with open(outfile, 'w') as f:
    for i,path in enumerate(paths):
        f.write('>Unitig%d\n'%(i))
        for j in range(len(path)-1):
            nodeA = path[j].lstrip("B")
            nodeB = path[j+1].lstrip("B")

            d =  g.get_edge_data(path[j],path[j+1])

            f.write('%s %s %s %s %d %d %d %d %d\n'%(nodeA.split('_')[0],nodeA.split('_')[1]  , nodeB.split('_')[0], 
                    nodeB.split('_')[1], -d['read_a_start_raw'] + d['read_a_end_raw'] - d['read_b_start_raw'] + d['read_b_end_raw'], 
                    d['read_a_start_raw'], d['read_a_end_raw'], d['read_b_start_raw'], d['read_b_end_raw']))



    f.close()




