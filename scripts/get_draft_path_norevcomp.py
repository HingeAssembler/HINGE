#!/usr/bin/env python

import sys
import os
from pbcore.io import FastaIO


def run(reader, writer):
    for i,record in enumerate(reader):
        if i%2 == 0:
            writer.writeRecord(record.header, record.sequence)


if __name__ == '__main__':
    iname, oname = sys.argv[1:3]
    reader = FastaIO.FastaReader(iname)
    writer = FastaIO.FastaWriter(oname)
    run(reader, writer)
