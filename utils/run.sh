#!/bin/bash

echo "Setting stuff up"
cur_fol=$PWD
cd ~/AwesomeAssembler && source utils/setup.sh
cd $cur_fol


echo "Running filter"
Reads_filter --las $1.las --db $1 --config ~/AwesomeAssembler/utils/nominal.ini

echo "Running hinging"

hinging --las $1.las --db $1 --config ~/AwesomeAssembler/utils/nominal.ini -o $1.edges

echo "Running Visualise"

python ~/AwesomeAssembler/scripts/Visualise_graph.py $1.edges.hinges hinge_list.txt

echo "Running Condense"

python ~/AwesomeAssembler/scripts/condense_graph.py $1.edges.hinges



echo "putting ground truth and condensing"