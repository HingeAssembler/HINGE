#!/bin/bash

echo "Setting stuff up"
cur_fol=$PWD
cd ~/AwesomeAssembler && source utils/setup.sh
cd $cur_fol


echo "Running filter"
Reads_filter --las $1.las --db $1 --config ~/AwesomeAssembler/utils/nominal.ini

echo "Running hinging"

hinging --las $1.las --db $1 --config ~/AwesomeAssembler/utils/nominal.ini -o $1.edges -x $1

echo "Running Visualise"

python ~/AwesomeAssembler/scripts/Visualise_graph.py $1.edges.hinges hinge_list.txt

echo "Running Condense"

python ~/AwesomeAssembler/scripts/condense_graph.py $1.edges.hinges

echo "Putting ground truth and condensing"
if [ -e "$1.mapping.1.json" ]
	then
	python ~/AwesomeAssembler/scripts/condense_graph_with_aln_json.py $1.edges.hinges $1.mapping.1.json
fi