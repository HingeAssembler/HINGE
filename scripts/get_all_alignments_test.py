import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact 
import interface_utils as util
import sys

import os
os.environ['PATH'] += ':/data/pacbio_assembly/AwesomeAssembler/DALIGNER'
#print os.popen("export").read()

path = '/data/pacbio_assembly/AwesomeAssembler/data/'
aln = []
aln = list(util.get_all_alignments2(path+'ecoli',path+'ecoli.las'))

print aln[0:5]