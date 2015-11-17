#!/usr/bin/python

import sys

def parse_qv(stream = sys.stdin):    
    with stream as f:
        for l in f:
            if l[0] == 'I':
                yield(l.split()[-1])
            else:
                pass
                
#for qv in parse_qv():
#    print qv
# do whatever you want to do with the reads
