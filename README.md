# AwesomeAssembler

![](https://magnum.travis-ci.com/Eureka22/AwesomeAssembler.svg?token=i41xfGcHb72GYFyZnvtg&branch=master)
## Design of the DB interface module

The DB interface module aims at providing fast and easy access to Dazzler DB and .las alignment file. The backend is a modification of LAshow module (in Fei's fork of DAligner), several filtering options will be added. The frondend is a python module for easy access and interactivity.

The python scripts take advantage of stdout to process large input and output. Functions are in a generator fashion.


## Usage

```bash
cd scripts
source setup.sh # setup the interface globally

cd <path of your database> 

./run_parse_read.py G 1 # assume your database is G.db
./run_parse_alignment.py G 1 

```

Sample output:

Reads:
```
('Sim/1/0_11053 RQ=0.850', 'cacagactcactccacactcgaatgtggatagggacactgtactcttcgatgcaaggaattggttaaccacgtttcgtggaaactcgggacgcttaggtaggctctcgtacagaccttcagcgaccacaagtcgattgaaagtgttctcatcaaggaacagactaagaaaaccccgacctgttactggaggggagtgaacagttagcgagtgccaccgagtcagaatcgagacgccctatctaaattgaagtaaatactgactaacggcgaagaccagctattgtagcaccccccgaccataactaaggggcgactacgcatggacggataaggcaccaacgtcgctaaccggtaaagtgcctcgcactggtcggacatcagatccttcgcctggctggtgatctaacaaaattatagcgaaaaggaccacgagactgataatggcttgctcgtattactactcgcagtgaaatcgtgccgtcctacggagtatactgttaaatctcacttggttcaccatgcgcggggcccctgactgcgaacggatcaactccagaaaagtaaggtattaacctagctcgttaatcggtgcaaacctgtcttaaactatgct...')
```


Alignments:
```
['n', 1, 32, 5121, 11053, 0, 5986, 1329, 60]
['c', 1, 116, 2982, 11053, 0, 8131, 1847, 82]
['c', 1, 119, 5806, 11053, 0, 5283, 1159, 53]
['c', 1, 193, 0, 1090, 8527, 9569, 230, 11]
['c', 1, 212, 1240, 10556, 0, 9326, 2091, 94]
['n', 1, 285, 6110, 11053, 0, 4994, 1118, 50]
['c', 1, 287, 0, 5351, 6510, 11790, 1182, 54]
['c', 1, 313, 0, 3897, 5130, 8978, 896, 39]
['c', 1, 354, 0, 4836, 7524, 12260, 1074, 49]
['n', 1, 491, 4493, 11053, 0, 6640, 1477, 67]
['n', 1, 492, 0, 3180, 7416, 10513, 689, 32]
['c', 1, 596, 859, 11053, 0, 10145, 2353, 103]
...
```



