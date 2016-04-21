#!/bin/bash
pwd=$PWD
cd $pwd/DAZZ_DB
make clean && make -j 8

cd $pwd/DALIGNER
make clean && make -j 8

cd $pwd/DASCRUBBER
make clean && make -j 8

cd $pwd
source setup.sh

mkdir $pwd/data
cd $pwd/data

#rm -rf *
rm G.*
simulator 1.0 -c50. >G.fasta
fasta2DB G G.fasta
DBsplit -s20 G
HPCdaligner G | csh -v 
rm G.*.G.*.las
LAmerge G G.*.las
DASqv -c50 G G.las

touch log.txt
LAInterface_test>log.txt
Consensus_test>log.txt
LAInterface_test_2DB>log.txt
