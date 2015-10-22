#!/usr/bin/python
import sys
import os
import subprocess
from parse_read import *
from parse_alignment import *

filename,filename2 = sys.argv[1:3]
alignmentname = sys.argv[3]
readarg = sys.argv[4]


stream = subprocess.Popen(["LA4Awesome", filename, filename2 , alignmentname ,readarg, '-F'],
                                  stdout=subprocess.PIPE, bufsize=1)

alignments = parse_alignment2(stream.stdout) # generator

d = {}
for alignment in alignments:
    if not d.has_key(alignment[2]):
        d[alignment[2]] = []
    d[alignment[2]].append([alignment[0],alignment[3],alignment[4], alignment[6], alignment[7]])
    
#print d

mapping = {}

for key,value in d.items():
    value.sort(key = lambda x:x[2]-x[1], reverse=True)
    aln = value[0]
    
    if aln[0] == 'n':
        mapping[str(key)] = (aln[1], aln[2])
        mapping[str(key)+'\''] = (aln[2], aln[1])
    else:
        mapping[str(key)] = (aln[2], aln[1])
        mapping[str(key)+'\''] = (aln[1], aln[2])
    
#print mapping
import ujson
ujson.dump(mapping,open('mapping.json','w'))
