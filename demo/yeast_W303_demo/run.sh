wget -nc -i https://gist.githubusercontent.com/pb-jchin/6359919/raw/9c172c7ff7cbc0193ce89e715215ce912f3f30e6/gistfile1.txt
dextract -o *.bax.h5


fasta2DB yeast m130605_000141_42207_c100515142550000001823076608221372_s1_p0.fasta 



DBsplit yeast

HPC.daligner yeast | bash -v

rm yeast.*.yeast.*.las
LAmerge yeast.las yeast.[0-9].las
DASqv -c100 yeast yeast.las

mkdir log


hinge filter --db yeast --las yeast --mlas -x yeast --config nominal.ini
hinge maximal --db yeast --las yeast --mlas -x yeast --config nominal.ini

hinge layout --db yeast --las yeast -x yeast --config nominal.ini -o yeast

hinge clip yeast.edges.hinges yeast.hinge.list demo

hinge draft-path $PWD yeast yeastdemo.G3.graphml
hinge draft --db yeast --las yeast.las --prefix yeast --config nominal.ini --out yeast.draft



hinge correct-head yeast.draft.fasta yeast.draft.pb.fasta draft_map.txt 
fasta2DB draft yeast.draft.pb.fasta

HPC.daligner yeast draft | bash -v 

rm draft.*.yeast.*.las
LAmerge draft.yeast.las draft.yeast.*.las

hinge consensus draft yeast draft.yeast.las yeast.consensus.fasta nominal.ini

hinge gfa $PWD yeast  yeast.consensus.fasta
