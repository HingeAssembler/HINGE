# AwesomeAssembler


## Design of the DB interface module

The DB interface module aims at providing fast and easy access to Dazzler DB and .las alignment file. The backend is a modification of LAshow module (in Fei's fork of DAligner), several filtering options will be added. The frondend is a python module for easy access and interactivity.

The python scripts take advantage of stdout to process large input and output. Functions are in a generator fashion.


## Usage

```
cd scripts
source setup.sh # setup the interface globally

cd <path of your database> 

./run_parse_read.py G 1 # assume your database is G.db
./run_parse_alignment.py G 1 

```

