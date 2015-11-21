# AwesomeAssembler

## Introduction 

AwesomeAssembler is an experimental long read assembler based on sparse string graph (|E|/|V| bounded). Now AwesomeAssembler is at research prototype stage.

## Pipeline Overview

AwesomeAssembler is an OLC(Overlap-Layout-Consensus) assembler. The idea of the pipeline is shown below. One significant difference from HGAP or Falcon pipeline is that it does not have a pre-assemble or read correction step. There are mainly two reasons that we don't want that step is that this step will sometimes collapse repeats, thus introduce errors to homologous regions. Secondly, this will throw away information, for example, half bases are thrown away (as for the ecoli dataset) in the preassembly step in falcon pipeline. 

![image](http://fxia.me/static/arc.png)

## Algorithm Details

### Reads filtering
Reads filtering filters reads that have long chimer in the middle, and short reads. 

### Layout 
For the layout step, currently there are two algorithms implemented. For both algorithms, one read and its reverse complement are kept as separated nodes (which needs to be changed later)

- `Greedy`. Implemented in `Greedy_best_ovl_bursty_resistant.cpp`. The Greedy works at all-approximate-repeats-bridged regime. For this algorithm, every best right extension for a read and its reverse complement are found. 
- `Not-so-Greedy`. Implemented in `NSG_v1.cpp`. NSG *aims* at working at all-triple-repeats are bridged regime. It is an augmented version for Greedy by adding a step to reduce false negative edges. 

### Pruning

Currently pruning is only available for `Greedy`. For the initial layout graph, it keeps removing degree-1 nodes, for a certain number of iterations. Implemented in `prune.py`. 

### Draft assembly

Draft assembly extracts bases from the backbone. First, it reduces the backbone to a linear structure by traversing the backbone assembly, as implemented in `draft_assembly.py`. Then, I use a unitig consensus algorithm to get the draft assembly from the linear structure, as in `Unitig_consensus.cpp`

### Consensus

After the draft assembly is obtained. Reads are mapped to the draft assembly and get the final assembly by doing a majority vote.

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
