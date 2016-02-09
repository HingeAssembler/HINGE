import sys 
import os



def run(inputfile,outputfile):


	nodes = {}
	arcs = {}

	with open(inputfile) as file:

	    for line in file:

			line_str = line[:-1]
			split_str = line_str.split(' ')

			node0 = int(split_str[0])
			node1 = int(split_str[1])

			# print node0,node1

			nodes[node0] = 1
			nodes[node1] = 1

			if node0 < node1:
				arcs[tuple([node0,node1])] = 1
			else:
				arcs[tuple([node1,node0])] = 1



	with open(outputfile, 'w') as fout:

		for node in nodes:
			fout.write("NODE "+str(node)+' 0 0 0 0 0\n')
			fout.write('AAA\n')
			fout.write('AAA\n')
			# print "NODE "+str(node)

		for arc in arcs:
			fout.write("ARC "+str(arc[0])+' '+str(arc[1])+' 0\n')

		



def main():

    run(sys.argv[1],sys.argv[2])
    return


if __name__ == '__main__':
    main()
