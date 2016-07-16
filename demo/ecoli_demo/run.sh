wget http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p4_filtered.fastq.gz
gunzip ecoli_p4_filtered.fastq.gz

seqtk seq -a ecoli_p4_filtered.fastq > reads.fasta
correct_head.py reads.fasta reads.pb.fasta map.txt
fasta2DB ecoli reads.pb.fasta

DBsplit ecoli

HPCdaligner ecoli | bash -v

rm ecoli.*.ecoli.*.las
LAmerge ecoli.las ecoli.+([[:digit:]]).las
DASqv -c100 ecoli ecoli.las

mkdir log

Reads_filter --db ecoli --las ecoli.las -x ecoli --config ~/HINGE/utils/nominal.ini
hinging --db ecoli --las ecoli.las -x ecoli --config ~/HINGE/utils/nominal.ini -o ecoli

pruning_and_clipping.py ecoli.edges.hinges ecoli.hinge.list demo

get_draft_path.py $PWD ecoli ecolidemo.G2.graphml
draft_assembly --db ecoli --las ecoli.las --prefix ecoli --config ~/HINGE/utils/nominal.ini --out ecoli.draft

correct_head.py ecoli.draft.fasta ecoli.draft.pb.fasta draft_map.txt 
fasta2DB draft ecoli.draft.pb.fasta
HPCmapper draft ecoli | bash -v 

rm draft.*.ecoli.*.las
LAmerge draft.ecoli.las draft.ecoli.*.las
consensus draft ecoli draft.ecoli.las ecoli.consensus.fasta ~/HINGE/utils/nominal.ini

get_consensus_gfa.py $PWD ecoli ecoldemo.G2.graphml ecoli.consensus.fasta
