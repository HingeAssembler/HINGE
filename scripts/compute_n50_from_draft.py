import sys
import os
import networkx as nx
from Bio import SeqIO



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


hinging_n50 = -1
hinging_comp_n50 = -1
hgap_n50 = -1

count = 0
count1 = 0
count2 = 0


data_dict = {}




fullpath = '/data/pacbio_assembly/pb_data/NCTC/'

for nctc_name in os.listdir(fullpath):
	if 'NCTC' not in nctc_name:
		continue

	mypath = fullpath+nctc_name

	if not os.path.isdir(mypath):
		continue

	mypath = mypath+'/'

	count += 1

	hinging_n50 = -1
	hinging_comp_n50 = -1
	hgap_n50 = -1


	data_dict[nctc_name] = []

	draft_file = [x for x in os.listdir(mypath) if 'draft.graphml' in x]

	try:

		# flname = sys.argv[1]
		g = nx.read_graphml(mypath+draft_file[0])

		contig_lengths = []
		component_lengths = []


		for u in g.nodes():
		    contig_lengths.append(len(g.node[u]['segment']))
		    
		for c in nx.weakly_connected_components(g):
		    # we use set() so that we cannot double-count a two reverse complementary contigs
		    # in the same component
		    component_lengths.append(sum(set([len(g.node[u]['segment']) for u in c])))
		    
		component_lengths = set(component_lengths)

		hinging_n50 = comp_n50(contig_lengths)
		hinging_comp_n50 = comp_n50(component_lengths)

		count1+=1

	except:
		pass

	# print "contig n50: "+str(comp_n50(contig_lengths))
	# print "component n50: "+str(comp_n50(component_lengths))


	hgap_file = [x for x in os.listdir(mypath) if 'hgap.fasta' in x]


	try:

		hgap_file = hgap_file[0]
		hgap_contigs = [len(x) for x in SeqIO.parse(open(mypath+hgap_file),'fasta')]

		hgap_n50 = comp_n50(hgap_contigs)

		count2+=1

	except:
		pass


	with open(mypath+nctc_name+'.n50','w') as f:
		f.write('hinging'+'\t'+str(hinging_n50)+'\n')
		f.write('hinging_comp'+'\t'+str(hinging_comp_n50)+'\n')
		f.write('hgap'+'\t'+str(hgap_n50)+'\n')

	data_dict[nctc_name] = [hinging_n50,hinging_comp_n50,hgap_n50]

	print count
	print count1
	print count2



with open(fullpath+'computed.n50','w') as f:
	for nctc_name in data_dict:
		vec = data_dict[nctc_name]
		f.write(nctc_name+'\t'+str(vec[0])+'\t'+str(vec[1])+'\t'+str(vec[2])+'\n')



