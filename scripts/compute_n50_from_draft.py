import sys
import networkx as nx



def comp_n50(contig_vec):
    if len(contig_vec) == 0:
        return 0
    sorted_lengths = sorted(contig_vec)

    total_length = sum(contig_vec)

    half_length = 0.5*total_length

    min_n50 = sorted_lengths[-1]
    max_n50 = 0

    for i in range(len(sorted_lengths)):
        sum_1 = sum(sorted_lengths[0:i+1])
        sum_2 = sum(sorted_lengths[i:])
        if sum_1 >= half_length and sum_2 >= half_length:
            min_n50 = min(sorted_lengths[i],min_n50)
            max_n50 = max(sorted_lengths[i],max_n50)

    return 0.5*(min_n50+max_n50)




flname = sys.argv[1]
g = nx.read_graphml(flname)


contig_lengths = []
component_lengths = []


for u in g.nodes():
    contig_lengths.append(len(g.node[u]['segment']))
    
for c in nx.weakly_connected_components(g):
    # we use set() so that we cannot double-count a two reverse complementary contigs
    # in the same component
    component_lengths.append(sum(set([len(g.node[u]['segment']) for u in c])))
    
component_lengths = set(component_lengths)



print "contig n50: "+str(comp_n50(contig_lengths))
print "component n50: "+str(comp_n50(component_lengths))



