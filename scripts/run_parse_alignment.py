#!/usr/bin/python
import sys
import os
import subprocess
from parse_read import *
from parse_alignment import *

filename = sys.argv[1]
readarg = sys.argv[2]


stream = subprocess.Popen(["LAshow", filename , filename ,readarg],
                                  stdout=subprocess.PIPE, bufsize=1)

alignments = parse_alignment(stream.stdout) # generator

for alignment in alignments:
    print alignment