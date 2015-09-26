#!/usr/bin/python

import sys
min_len_aln = 1000


with sys.stdin as f:
    for l in f:
        l = l.strip().split()
        if len(l) != 2:
            continue

        read_id = l[0]
        seq = l[1]
        
        print read_id,seq
        
        #if len(seq) > max_len:
        #    seq = seq[:max_len-1]
        
        if read_id not in ("+", "-", "*"):
            if len(seq) >= min_len_aln:
                if len(seqs) == 0:
                    seqs.append(seq) #the "seed"
                    seed_id = l[0]
                if read_id not in read_ids: #avoidng using the same read twice. seed is used again here by design
                    seqs.append(seq)
                    read_ids.add(read_id)
        elif l[0] == "+":
            if len(seqs) >= min_cov_aln:
                seqs = seqs[:1] + sorted(seqs[1:], key=lambda x: -len(x))
                yield (seqs[:max_n_read], seed_id, config) 
            #seqs_data.append( (seqs, seed_id) ) 
            seqs = []
            read_ids = set()
            seed_id = None
        elif l[0] == "*":
            seqs = []
            read_ids = set()
            seed_id = None
        elif l[0] == "-":
            #yield (seqs, seed_id)
            #seqs_data.append( (seqs, seed_id) )
            break
