#!/usr/bin/python

import sys

def parse_read(stream = sys.stdin):    
    sid = ''
    seq = ''
    with stream as f:
        for l in f:
            if l[0] == '>':
                if sid == '':
                    sid = l[1:].strip()
                else:
                    tsid = sid
                    tseq = seq
                    seq = ''
                    sid = l[1:].strip()
                    yield (tsid,tseq)
                    
            else:
                seq += l.strip()
            
        yield(sid,seq)
                
                
#for read in parse_read():
#    print read
# do whatever you want to do with the reads
