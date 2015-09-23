#!/usr/bin/python

import sys
import os
import subprocess
from parse_read import *

filename = sys.argv[1]
readarg = sys.argv[2]


stream = subprocess.Popen(["DBshow", filename ,readarg],
                                  stdout=subprocess.PIPE,bufsize=1)

reads = parse_read(stream.stdout) # generator

for read in reads:
    print read

#print result