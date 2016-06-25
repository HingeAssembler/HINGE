#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ipywidgets.widgets import interact
import interface_utils as util
import sys

import os
os.environ['PATH'] += ':~/AwesomeAssembler/DALIGNER'
#print os.popen("export").read()

left = int(sys.argv[1])
right = int(sys.argv[2])

ref = sys.argv[3]
read = sys.argv[4]
las = sys.argv[5]
contig = sys.argv[6]
length_th = int(sys.argv[7])

#path = '/data/pacbio_assembly/AwesomeAssembler/data/ecoli/'
aln = []

#bb = []
#with open('ecoli.linear.edges') as f:
#    for line in f:
#        e = line.split(" ")[0]
#        if e[-1] == '\'':
#            e = e[:-1]
#
#        bb.append(int(e))
#
#print bb
#
#bb = set(bb)


for i,item in enumerate(util.get_alignments_mapping3(ref, read, las, contig)):
    if i%2000 == 0:
        print i, item

    if item[3] >= left and item[4] <= right and item[4] - item[3] > length_th:
        aln.append(item)


print 'number:',len(aln)
aln.sort(key = lambda x:x[2])

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

alns.sort(key = lambda x:min([item[3] for item in x]))


plt.figure(figsize = (15,10))
plt.axes()
#plt.gca().axes.get_yaxis().set_visible(False)
#l = aln[0][5]
tip = (right-left)/5000
ed = (right-left)/2000
grid_size = 1.0
plt.xlim(left-2000,right+2000)
plt.ylim(-5,num*grid_size)

points = [[left,0], [right,0], [right+tip,grid_size/4], [right,grid_size/2], [left,grid_size/2]]
#rectangle = plt.Rectangle((0, 0), l, 5, fc='r',ec = 'none')
polygon = plt.Polygon(points,fc = 'r', ec = 'none', alpha = 0.6)
plt.gca().add_patch(polygon)

dotted_line = plt.Line2D((left, left), (0, num*grid_size ),ls='-.')
plt.gca().add_line(dotted_line)

dotted_line2 = plt.Line2D((right, right), (0, num*grid_size ),ls='-.')
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

        #if item[2] in bb:
        #    polygon = plt.Polygon(points,fc = 'r', ec = 'none', alpha = 0.8)
        #else:
        #    polygon = plt.Polygon(points,fc = 'b', ec = 'none', alpha = 0.6)

        polygon.set_url("http://shannon.stanford.edu:5000/aln" + str(item[2]+1) + ".pdf")
        plt.gca().add_patch(polygon)

        if points_end != []:
            polygon2 = plt.Polygon(points_end,fc = 'g', ec = 'none', alpha = 0.6)
            plt.gca().add_patch(polygon2)

        if points_start != []:
            polygon2 = plt.Polygon(points_start,fc = 'g', ec = 'none', alpha = 0.6)
            plt.gca().add_patch(polygon2)

plt.savefig('mapping/map.' + str(left) +'_'+ str(right)+ '.svg')
