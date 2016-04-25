#!/usr/bin/python

import sys
import os


bact_name = "ecoli"

if len(sys.argv) >= 2:
	bact_name = sys.argv[1]

st_point = 0
if len(sys.argv) >= 3:
	st_point = int(sys.argv[2])




dextract_command = 'dextract -o '
h5list = []

fasta_name = ''

for filename in os.listdir('.'):
	if len(filename.split('.')) > 2 and filename.split('.')[-2] == 'bax' and filename.split('.')[-1] == 'h5':
		h5list.append(filename)
		dextract_command = dextract_command + ' ' + filename
		if filename.split('.')[-3] == '1':
			fasta_name = filename.split('.1.bax.h5')[0] + '.fasta'


if st_point <= 0:
	print dextract_command
	os.system(dextract_command)

if st_point <= 1:
	print "fasta2DB "+bact_name+' '+fasta_name
	os.system("fasta2DB "+bact_name+' '+fasta_name)

if st_point <= 2:
	print "DBsplit -x500 -s100 "+bact_name
	os.system("DBsplit -x500 -s100 "+bact_name)

if st_point <= 3:
	print "HPCdaligner -dal4 -t16 -e.7 -l500 -s100 "+bact_name+" | zsh"
	os.system("HPCdaligner -dal4 -t16 -e.7 -l500 -s100 "+bact_name+" | zsh")

if st_point <= 4:
	print "rm "+bact_name+".*."+bact_name+".*"
	os.system("rm "+bact_name+".*."+bact_name+".*")

if st_point <= 5:
	print "LAmerge "+bact_name+".las "+bact_name+".*.las"
	os.system("LAmerge "+bact_name+".las "+bact_name+".*.las")

if st_point <= 6:
	print "rm "+bact_name+".*.las"
	os.system("rm "+bact_name+".*.las")

if st_point <= 7:
	os.system("mkdir log")

if st_point <= 8:
	print "DASqv -c100 "+bact_name+" "+bact_name+".las"
	os.system("DASqv -c100 "+bact_name+" "+bact_name+".las")

if st_point <= 9:
	print "Reads_filter --db "+bact_name+" --las "+bact_name+".las -x "+bact_name+" --config ~/AwesomeAssembler/utils/nominal.ini"
	os.system("Reads_filter --db "+bact_name+" --las "+bact_name+".las -x "+bact_name+" --config ~/AwesomeAssembler/utils/nominal.ini")

if st_point <= 10:
	print "hinging --db "+bact_name+" --las "+bact_name+".las -x "+bact_name+" --config ~/AwesomeAssembler/utils/nominal.ini -o "+bact_name
	os.system("hinging --db "+bact_name+" --las "+bact_name+".las -x "+bact_name+" --config ~/AwesomeAssembler/utils/nominal.ini -o "+bact_name)

if st_point <= 11:
	print "python ~/AwesomeAssembler/scripts/pruning_and_clipping.py "+bact_name+".edges.hinges "+bact_name+".hinge.list A"
	os.system("python ~/AwesomeAssembler/scripts/pruning_and_clipping.py "+bact_name+".edges.hinges "+bact_name+".hinge.list A")


