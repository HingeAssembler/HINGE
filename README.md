# AwesomeAssembler

CI Status: ![image](https://magnum.travis-ci.com/Eureka22/AwesomeAssembler.svg?token=i41xfGcHb72GYFyZnvtg&branch=master)

## Introduction 

AwesomeAssembler is an experimental long read assembler based on sparse string graph (|E|/|V| bounded). Now AwesomeAssembler is at research prototype stage.

## Pipeline Overview

AwesomeAssembler is an OLC(Overlap-Layout-Consensus) assembler. The idea of the pipeline is shown below. One significant difference from HGAP or Falcon pipeline is that it does not have a pre-assembly or read correction step. There are mainly two reasons that we don't want that step. The first is that this step will sometimes collapse repeats, thus introduce errors to homologous regions. Secondly, this will throw away information, for example, half bases are thrown away (as for the ecoli dataset) in the preassembly step in falcon pipeline. 

![image](http://fxia.me/static/arc.png)

## Algorithm Details

### Reads filtering
Reads filtering filters reads that have long chimer in the middle, and short reads. 

### Layout 
For the layout step, currently there are two algorithms implemented. For both algorithms, one read and its reverse complement are kept as separated nodes (which needs to be changed later)

- `Greedy`. Implemented in `Greedy_best_ovl_bursty_resistant.cpp`. The Greedy works at all-approximate-repeats-bridged regime. For this algorithm, every best right extension for a read and its reverse complement are found. 
- `Not-so-Greedy`. Implemented in `NSG_v1.cpp`. NSG *aims* at working at all-triple-repeats are bridged regime. It is an augmented version for Greedy by adding a step to reduce false negative edges. 

### Pruning

Currently pruning is only available for `Greedy`. For the initial layout graph, it keeps removing degree-1 nodes, for a certain number of iterations. Implemented in `prune.py`. One output is a graphml file which is the graph representation of the backbone. It can be analyzed and visualized, etc. 

### Draft assembly

Draft assembly extracts bases from the backbone. First, it reduces the backbone to a linear structure by traversing the backbone assembly, as implemented in `draft_assembly.py`. Then, I use a unitig consensus algorithm to get the draft assembly from the linear structure, as in `Unitig_consensus.cpp`

### Consensus

After the draft assembly is obtained. Reads are mapped to the draft assembly and get the final assembly by doing a majority vote. In `Draft_consensus.cpp`.

## Parameters

In the pipeline described above, most programs not only takes the input file and output file as arguments, but also require a configuration file in ini format. This consists parameters for each step in the pipeline, and their usage and effects are explained below:


###[filter]
- length_threshold = 6500; // Length threshold for reads to be considered in the backbone
- quality_threshold = 0.23; // Quality threshold for edges to be considered in the backbone 
- n_iter = 2; // iterations of filtering, the filtering needs several iterations, because when filter reads, you got rid of some edges; when filter edges, you got rid of some reads (if the last edge is filtered.) Typically 2-3 iterations will be enough.
- aln_threshold = 2500; // Length of alignment for edges to be considered in the backbone
- min_cov = 5; // Minimal coverage for a segment to be considered not chimer/adaptor
- cut_off = 200; // A parameter for identifying long chimer in the middle of a read
- theta = 300; // A parameter for tolerance of the overhang length when looking for right extension.


###[running]
- n_proc = 12; // number of CPUs for layout step

###[draft]
- min_cov = 10; //obsolete
- trim = 200; //obsolete
- edge_safe = 100; //obsolete
- tspace = 900; //space between new "trace points"


###[consensus]
- min_length = 2000; // Minimal length of reads used for final consensus
- trim_end = 200; // Trim ends for alignments for final consensus
- best_n = 1; // If one read has multiple alignments with the bacbone assembly, choose the longest n segments for consensus.
- quality_threshold = 0.23; // alignment quality threshold

# Installation

This software is still at prototype stage so it is not well packaged, however it is designed in a modular flavor so different combinations of methods can be tested. 

Installing the software is very easy. 

```
git clone https://github.com/Eureka22/AwesomeAssembler.git
git submodule init
git submodule update
./build.sh
```

# Running

In order to call the programs from anywhere, I suggest one export the directory of binary file to system environment, you can do that by using the script `setup.sh`.

A demo run for assembling the ecoli genome is the following:

```
source setup.sh %I've changed the setup.sh, as it gave errors while running. -GK
mkdir data/ecoli
cd data/ecoli
# reads.fasta should be in data/ecoli
fasta2DB ecoli reads.fasta
DBsplit -x500 -s100 ecoli     
HPCdaligner -dal4 -t16 -e.7 -l500 -s100 ecoli | zsh
# alternatively, you can put output of HPCdaligner to a bash file and edit it to support 
rm ecoli.*.ecoli.*
LAmerge ecoli.las ecoli.*.las
rm ecoli.*.las # we only need ecoli.las

Greedy_best_ovl ecoli ecoli.las ecoli_greedy.edges greedy.ini
prune.py ecoli_greedy.edges
draft_assembly.py ecoli.edges
Unitig_consensus ecoli ecoli.las ecoli.linear.edges ecoli.draft.fasta greedy.ini

correct_head.py ecoli.draft.fasta ecoli.draft.pb.fasta 
fasta2DB draft ecoli.draft.pb.fasta
HPCmapper -e.73 draft ecoli | zsh -v 
LAmerge draft.ecoli.las draft.ecoli.*.las
rm draft.ecoli.*.las
Draft_consensus draft ecoli draft.ecoli.las ecoli.consensus.fasta greedy.ini 
# final consensus is in ecoli.consensus.fasta
```

## Debugging

### showing ground truth on graph
Some programs are for debugging and oberservation. For example, one can get the ground truth by mapping reads to reference and get `ecoli.ecoli.ref.las`.

This `las` file can be parsed to json file for other programs to use. 

```
run_mapping.py ecoli ecoli.ref ecoli.ecoli.ref.las 1-$ 
```

In the prune step, if `ecoli.mapping.json` exists, the output `graphml` file will contain the information of ground truth. 

### drawing alignment graphs and mapping graphs
Draw a read, for example 60947, and output figure to `sample` folder (need plus 1 as LAshow counts from 1):

```
draw2.py ecoli ecoli.las 60948 sample 100
```

Draw pileup on draft assembly, given a region(start,end):

```
draw2_pileup_region.py  3600000 4500000 
```

# Preliminary results:
![image](http://fxia.me/assets/img/awe.png)
For ecoli 160X dataset, finished assembly can be achieved. 


# Limitations

- Only tested on high coverage and microbe datasets
- Keeping one read and its reverse complement as two nodes lose a lot of information
- More work needed to make NSG and Z-Pruning working. 
