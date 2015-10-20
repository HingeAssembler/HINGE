#!/usr/bin/python
import sys, os
from pbcore.io import FastaIO

def run(reader, writer):
    for record in reader:
        seq_length = len(record.sequence)
        zmw = 0
        #bounds = record.header.split('/')[-1]
        #start, end = [int(k) for k in bounds.split('_')]
        start = 0
        new_end = start + seq_length

        new_header = "m000_000/{zmw}/{start}_{end}".format(zmw=zmw, start=start, end=new_end)

        writer.writeRecord(new_header, record.sequence)

def main(iname, ofile):
    reader = FastaIO.FastaReader(iname)
    writer = FastaIO.FastaWriter(ofile)
    run(reader, writer)

if __name__ == '__main__':
    iname, oname = sys.argv[1:3]
    ofile = open(oname, 'w')
    try:
        main(iname, ofile)
    except:
        # clean up (for make)
        ofile.close()
        os.unlink(oname)
        raise