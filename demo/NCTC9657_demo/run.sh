hinge correct-head NCTC9657_reads.fasta reads.pb.fasta map.txt
fasta2DB NCTC9657 reads.pb.fasta

DBsplit NCTC9657

HPC.daligner NCTC9657 | bash -v

rm NCTC9657.*.NCTC9657.*.las
LAmerge NCTC9657.las NCTC9657.[0-9].las
DASqv -c100 NCTC9657 NCTC9657.las

mkdir log



hinge filter --db NCTC9657 --las "NCTC9657.*.las" -x NCTC9657 --config ../../utils/nominal.ini
hinge layout --db NCTC9657 --las NCTC9657.las -x NCTC9657 --config ../../utils/nominal.ini -o NCTC9657

hinge clip NCTC9657.edges.hinges NCTC9657.hinge.list demo

hinge draft-path $PWD NCTC9657 NCTC9657demo.G2.graphml
hinge draft --db NCTC9657 --las NCTC9657.las --prefix NCTC9657 --config ../../utils/nominal.ini --out NCTC9657.draft






hinge correct-head NCTC9657.draft.fasta NCTC9657.draft.pb.fasta draft_map.txt 
fasta2DB draft NCTC9657.draft.pb.fasta
HPC.daligner NCTC9657 draft | bash -v 

# rm draft.*.NCTC9657.*.las
# LAmerge draft.NCTC9657.las draft.NCTC9657.*.las

hinge consensus draft NCTC9657 draft.NCTC9657.las NCTC9657.consensus.fasta ../../utils/nominal.ini

hinge gfa $PWD NCTC9657 NCTC9657.consensus.fasta





