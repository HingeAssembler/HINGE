#!/usr/bin/python
import sys, os
from pbcore.io import FastaIO

def run(reader, writer, lookupfile):
    with open (lookupfile,'w') as f:
        for i,record in enumerate(reader):
            seq_length = len(record.sequence)
            zmw = i+1
            old_header=record.header
            #bounds = record.header.split('/')[-1]
            #start, end = [int(k) for k in bounds.split('_')]
            start = 0
            new_end = start + seq_length

            new_header = "m000_000/{zmw}/{start}_{end}".format(zmw=zmw, start=start, end=new_end)
            f.write(old_header+'\t'+new_header+'\n')

            writer.writeRecord(new_header, record.sequence)

def main(iname, ofile, lookupfile):
    reader = FastaIO.FastaReader(iname)
    writer = FastaIO.FastaWriter(ofile)
    run(reader, writer,lookupfile)

if __name__ == '__main__':
    iname, oname, lookupfile = sys.argv[1:4]
    ofile = open(oname, 'w')
    try:
        main(iname, ofile, lookupfile)
    except:
        # clean up (for make)
        ofile.close()
        os.unlink(oname)
        raise
