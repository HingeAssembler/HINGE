#!/usr/bin/python
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

n = (sys.argv[1])
rst = []
with open(n) as f:
    for line in f:
        tmp = line.strip().split()
        t1 = tmp[0]
        if t1[-1] == '\'':
            t1 = t1[:-1]
        t2 = tmp[1]
        if t2[-1] == '\'':
            t2 = t2[:-1]
        rst.append(int(t1)+1)
        #rst.append(int(t2)+1)

rst = range(1,1399)
path = '/data/pacbio_assembly/AwesomeAssembler/data/'
aln = []
for i,n in enumerate(rst):
    print i
    for item in util.get_alignments_mapping(path+'corr', path + 'ecoli.ref', path+'corr.ecoli.ref.las', [n]):
        aln.append(item)

print aln[0:20]

#aln.sort(key = lambda x:x[2])

alns = []
current_b = aln[0][2]
aln_group = []

for item in aln:
    if current_b != item[2]:
        alns.append(aln_group)
        aln_group = []
        aln_group.append(item)
        current_b = item[2]
    else:
        aln_group.append(item)

num = len(alns)
print len(aln), len(alns)

#print [len(item) for item in alns]
#print [item[0:3] for item in aln]

#alns.sort(key = lambda x:min([item[3] for item in x]))

#size_chunk = num/grid_size
#for i in range(grid_size):
#    aln[i*size_chunk:min((i+1)*size_chunk, num)] = sorted(aln[i*size_chunk:min((i+1)*size_chunk, num)],key = lambda x: x[4]-x[3] ,reverse=True)

plt.figure(figsize = (15,10))
plt.axes()
#plt.gca().axes.get_yaxis().set_visible(False)
l = aln[0][5]
tip = l/5000
ed = l/2000
grid_size = 1.0
plt.xlim(-2000,l+2000)
plt.ylim(-5,num*grid_size)

points = [[0,0], [l,0], [l+tip,grid_size/4], [l,grid_size/2], [0,grid_size/2]]
#rectangle = plt.Rectangle((0, 0), l, 5, fc='r',ec = 'none')
polygon = plt.Polygon(points,fc = 'r', ec = 'none', alpha = 0.6)
plt.gca().add_patch(polygon)

dotted_line = plt.Line2D((0, 0), (0, num*grid_size ),ls='-.')
plt.gca().add_line(dotted_line)

dotted_line2 = plt.Line2D((l, l), (0, num*grid_size ),ls='-.')
plt.gca().add_line(dotted_line2)

for i,aln_group in enumerate(alns):
    for item in aln_group:
        abpos = item[3]
        aepos = item[4]
        bbpos = item[6]
        bepos = item[7]
        blen = item[8]
        strand = item[0]
        points_start = []
        points_end = []

        if strand == 'n':
            points = [[abpos, (i+1)*grid_size], [aepos, (i+1)*grid_size], [aepos + tip, (i+1)*grid_size + grid_size/4], [aepos, (i+1)*grid_size+grid_size/2], [abpos, (i+1)*grid_size+grid_size/2]]
            if (bepos < blen):
                points_end = [[aepos, (i+1)*grid_size], [aepos + tip, (i+1)*grid_size + grid_size/4], [aepos, (i+1)*grid_size+grid_size/2], [aepos+ed, (i+1)*grid_size+grid_size/2], [aepos + ed+ tip, (i+1)*grid_size + grid_size/4],  [aepos+ed, (i+1)*grid_size]]
            if (bbpos > 0):
                points_start = [[abpos, (i+1)*grid_size], [abpos, (i+1)*grid_size+grid_size/2], [abpos-ed, (i+1)*grid_size+grid_size/2], [abpos-ed, (i+1)*grid_size]]
        else:
            points = [[abpos, (i+1)*grid_size], [aepos, (i+1)*grid_size], [aepos, (i+1)*grid_size+grid_size/2], [abpos, (i+1)*grid_size+grid_size/2], [abpos - tip, (i+1)*grid_size + grid_size/4]]
            if (bepos < blen):
                points_end = [[aepos, (i+1)*grid_size],  [aepos, (i+1)*grid_size+grid_size/2], [aepos+ed, (i+1)*grid_size+grid_size/2], [aepos+ed, (i+1)*grid_size]]
            if (bbpos > 0):
                points_start = [[abpos, (i+1)*grid_size],[abpos-tip, (i+1)*grid_size+grid_size/4], [abpos, (i+1)*grid_size+grid_size/2], [abpos-ed, (i+1)*grid_size+grid_size/2],[abpos-ed-tip, (i+1)*grid_size+grid_size/4], [abpos-ed, (i+1)*grid_size]]

        polygon = plt.Polygon(points,fc = 'b', ec = 'none', alpha = 0.6)
        polygon.set_url("http://shannon.stanford.edu:5000/aln" + str(item[2]+1) + ".pdf")
        plt.gca().add_patch(polygon)

        if points_end != []:
            polygon2 = plt.Polygon(points_end,fc = 'g', ec = 'none', alpha = 0.6)
            plt.gca().add_patch(polygon2)

        if points_start != []:
            polygon2 = plt.Polygon(points_start,fc = 'g', ec = 'none', alpha = 0.6)
            plt.gca().add_patch(polygon2)

plt.savefig('mapping/map.' + str(n)+ '.svg')
