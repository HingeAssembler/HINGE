#wget http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p4_filtered.fastq.gz
#gunzip ecoli_p4_filtered.fastq.gz

#seqtk seq -a ecoli_p4_filtered.fastq > reads.fasta
hinge correct-head reads.fasta reads.pb.fasta map.txt
fasta2DB ecoli reads.pb.fasta


DBsplit ecoli

HPC.daligner ecoli | bash -v

rm ecoli.*.ecoli.*.las
LAmerge ecoli.las ecoli.[0-9].las
DASqv -c100 ecoli ecoli.las

mkdir log



hinge filter --db ecoli --las "ecoli.*.las" -x ecoli --config ../../utils/nominal.ini
hinge layout --db ecoli --las ecoli.las -x ecoli --config ../../utils/nominal.ini -o ecoli

hinge clip ecoli.edges.hinges ecoli.hinge.list demo

hinge draft-path $PWD ecoli ecolidemo.G2.graphml
hinge draft --db ecoli --las ecoli.las --prefix ecoli --config ../../utils/nominal.ini --out ecoli.draft



hinge correct-head ecoli.draft.fasta ecoli.draft.pb.fasta draft_map.txt 
fasta2DB draft ecoli.draft.pb.fasta

HPC.daligner ecoli draft | bash -v 

#rm draft.*.ecoli.*.las
#LAmerge draft.ecoli.las draft.ecoli.*.las

hinge consensus draft ecoli draft.ecoli.las ecoli.consensus.fasta ../../utils/nominal.ini

hinge gfa $PWD ecoli ecolidemo.G2.graphml ecoli.consensus.fasta
