wget -nv http://files.pacb.com/datasets/secondary-analysis/ecoli-k12-P4C2-20KSS/ecoliK12.tar.gz
tar -zxf ecoliK12.tar.gz

dextract -o ecoliK12/Analysis_Results/*.bax.h5
fasta2DB ecoli m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.fasta


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

hinge gfa $PWD ecoli  ecoli.consensus.fasta
