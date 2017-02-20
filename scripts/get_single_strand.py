#!/usr/bin/env python

#usage python get_single_strand.py <in-fasta> <out-fasta>

from pbcore.io import FastaIO
import sys

flpath = sys.argv[1]
outpath = sys.argv[2]
writer = FastaIO.FastaWriter(outpath)
reader = FastaIO.FastaReader(flpath)
j = 0
for i,record in enumerate(reader):
    if j%2 == 0:
        writer.writeRecord('Consensus'+str(j), record.sequence)
        j+=1