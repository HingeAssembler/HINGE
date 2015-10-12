#!/bin/bash
pwd=$PWD
cd $pwd/DAZZ_DB
make clean && make -j 8

cd $pwd/DALIGNER
make clean && make -j 8

cd $pwd
source setup.sh

mkdir $pwd/data
cd $pwd/data

#rm -rf *
simulator 1.0 -c50. >G.fasta
fasta2DB G G.fasta 
HPCdaligner -mdust -t5 G | csh -v 

LAInterface_test
