import sys

chr_lengths={}

ground_truth=sys.argv[1]
graph_file=sys.argv[2]
out_file=sys.argv[2]+'.clipped'

with open(ground_truth) as f:
    for line in f:
        m = map(int, line.strip().split())
        chr_lengths.setdefault(m[1],0)
        chr_lengths[m[1]]= max(chr_lengths[m[1]], max(m[2],m[3]))

CHR_THR=20000

reads_to_kill=set()

with open(ground_truth) as f:
    for line in f:
        m = map(int, line.strip().split())
        read_left=min(m[2],m[3])
        read_right=max(m[2],m[3])
        read_chr=m[1]
        if read_left < CHR_THR:
            reads_to_kill.add(m[0])
        if read_right > chr_lengths[read_chr] - CHR_THR:
            reads_to_kill.add(m[0])

with open(graph_file) as f:
    with open(out_file, 'w') as g:
        for line in f:
            line1=line.split()
            if int(line1[0])in reads_to_kill or int(line1[1]) in reads_to_kill:
                continue
            g.write(line)
