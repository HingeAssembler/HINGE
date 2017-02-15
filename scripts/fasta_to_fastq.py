#!/usr/bin/env python
"""
Convert FASTA to FASTQ file with a static

Usage:
$ ./fasta_to_fastq NAME.fasta NAME.fastq
"""

import sys, os
from Bio import SeqIO

# Get inputs
fa_path = sys.argv[1]
fq_path = sys.argv[2]

# make fastq
with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
    for record in SeqIO.parse(fasta, "fasta"):
        record.letter_annotations["phred_quality"] = [40] * len(record)
        SeqIO.write(record, fastq, "fastq")
