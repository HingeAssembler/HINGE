#!/usr/bin/python
import sys
import os
import subprocess
from parse_read import *
from parse_alignment import *

filename = sys.argv[1]
alignmentname = sys.argv[2]
readarg = sys.argv[3]


stream = subprocess.Popen(["LA4Awesome", filename, filename , alignmentname ,readarg],
                                  stdout=subprocess.PIPE, bufsize=1)

alignments = parse_alignment2(stream.stdout) # generator

for alignment in alignments:
    print alignment