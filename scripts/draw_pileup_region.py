#!/usr/bin/env python
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



covy = np.zeros((right - left, ))
for item in aln:
    covy[item[3] - left : item[4] - left] += 1

covx = np.arange(left, right)


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

#num = len(alns)
num = len(aln)

print len(aln), len(alns)

alns.sort(key = lambda x:min([item[3] for item in x]))



fig = plt.figure(figsize = (15,10))
plt.axes()
ax1 = plt.subplot2grid((6,6), (0, 0), colspan=6, rowspan=4)
ax2 = plt.subplot2grid((6,6), (4, 0), colspan=6, rowspan=1, sharex = ax1)


#plt.gca().axes.get_yaxis().set_visible(False)
#l = aln[0][5]
tip = (right-left)/5000
ed = (right-left)/2000
grid_size = 1.0
ax1.set_xlim(left-2000,right+2000)
ax1.set_ylim(-5,num*grid_size)

points = [[left,0], [right,0], [right+tip,grid_size/4], [right,grid_size/2], [left,grid_size/2]]
#rectangle = plt.Rectangle((0, 0), l, 5, fc='r',ec = 'none')
polygon = plt.Polygon(points,fc = 'r', ec = 'none', alpha = 0.6)
ax1.add_patch(polygon)

dotted_line = plt.Line2D((left, left), (0, num*grid_size ),ls='-.')
ax1.add_line(dotted_line)

dotted_line2 = plt.Line2D((right, right), (0, num*grid_size ),ls='-.')
ax1.add_line(dotted_line2)

alns_all = []
for item in alns:
    for aln in item:
        alns_all.append([aln])

alns_all.sort(key = lambda x:min([item[3] for item in x]))


for i,aln_group in enumerate(alns_all):
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
        polygon = plt.Polygon(points,fc = 'b', ec = 'none', alpha = 0.6)

        polygon.set_url("http://shannon.stanford.edu:5000/aln" + str(item[2]+1) + ".pdf")
        ax1.add_patch(polygon)

        if points_end != []:
            polygon2 = plt.Polygon(points_end,fc = 'g', ec = 'none', alpha = 0.6)
            ax1.add_patch(polygon2)

        if points_start != []:
            polygon2 = plt.Polygon(points_start,fc = 'g', ec = 'none', alpha = 0.6)
            ax1.add_patch(polygon2)


ax2.plot(covx, covy)
plt.xlabel('position')
ax1.set_ylabel('pile-o-gram')
ax2.set_ylabel('coverage')


plt.savefig('mapping/map.' + str(contig) + '_' + str(left) +'_'+ str(right)+ '.svg')
