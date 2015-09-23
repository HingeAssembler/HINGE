#!/usr/bin/python

import sys
import os
import subprocess
from parse_read import *
from parse_alignment import *

#filename = sys.argv[1]
#readarg = sys.argv[2]

def get_reads(filename, readlist):
    stream = subprocess.Popen(["DBshow", filename , ' '.join(map(str,readlist))],
                                      stdout=subprocess.PIPE,bufsize=1)
    reads = parse_read(stream.stdout) # generator
    return reads

def get_alignments(filename, readlist):
    stream = subprocess.Popen(["LAshow", filename,filename , ' '.join(map(str,readlist))],
                                      stdout=subprocess.PIPE,bufsize=1)
    alignments = parse_alignment(stream.stdout) # generator
    return alignments

def get_all_reads(filename):
    stream = subprocess.Popen(["DBshow", filename],
                                      stdout=subprocess.PIPE,bufsize=1)
    reads = parse_read(stream.stdout) # generator
    return reads

def get_all_alignments(filename):
    stream = subprocess.Popen(["LAshow", filename, filename ],
                                      stdout=subprocess.PIPE,bufsize=1)
    alignments = parse_alignment(stream.stdout) # generator
    return alignments


# test   
#for item in get_reads('G',[1]):
#    print item
    
#for item in get_alignments('G',[1]):
#    print item

#for item in get_all_alignments('G'):
#    print item