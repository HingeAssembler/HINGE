% HINGE(1)
%
% October 2016

# NAME

hinge - assembler for long-read sequencing data

# SYNOPSIS

**hinge** {**subcommand**} *options* *files*

# OPTIONS

Subcommands are described below.
Run each subcommand without arguments for usage information.

**filter**
:    filter out short reads and long chimeric reads.

**maximal**
:    get maximal reads.

**layout**
:    generate a layout for assembly

**clip**
:    prune and clip output of the **layout** command

**draft-path**
:    get assembly graph as list of nodes

**draft**
:    construct draft assembly

**correct-head**
:    convert fasta file to daligner-specific format

**consensus**
:    construct consensus sequence

**gfa**
:    Create a graphical fragment assembly file from **consensus** output

**visualize** **visualise**

