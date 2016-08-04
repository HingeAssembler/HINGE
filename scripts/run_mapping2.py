#!/usr/bin/env python
import sys
import os
import subprocess
from parse_read import *
from parse_alignment import *

filename,filename2 = sys.argv[1:3]
alignmentname = sys.argv[3]
readarg = sys.argv[4]
k = int(sys.argv[5])


stream = subprocess.Popen(["LA4Awesome", filename, filename2 , alignmentname ,readarg],
                                  stdout=subprocess.PIPE, bufsize=1)

alignments = parse_alignment2(stream.stdout) # generator

d = {}
for alignment in alignments:
    if not d.has_key(alignment[2]):
        d[alignment[2]] = []
    d[alignment[2]].append([alignment[0],alignment[3],alignment[4], alignment[6], alignment[7], alignment[1]])
    
#print d

mapping = {}

for key,value in d.items():
    value.sort(key = lambda x:x[2]-x[1], reverse=True)
    alns = value[:k]
    max_val=alns[0][2]-alns[0][1]
    for aln in alns:
        if aln[2]-aln[1] > max_val/2.:
            if not mapping.has_key(str(key)):
                mapping[str(key)] = [(aln[1], aln[2],aln[-1], 1-int(aln[0] == 'n'))]
                # mapping[str(key)+'\''] = [(aln[2], aln[1],aln[-1], int(aln[0] == 'n'))]
            else:
                mapping[str(key)].append((aln[1], aln[2],aln[-1], 1-int(aln[0] == 'n')))
                # mapping[str(key)+'\''].append((aln[2], aln[1],aln[-1], int(aln[0] == 'n')))
            
#print mapping
import ujson
ujson.dump(mapping,open(filename2+'.mapping.'+str(k)+'.json','w'))
