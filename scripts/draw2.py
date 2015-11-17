#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from ipywidgets.widgets import interact 
import interface_utils as util
import sys
import os
import linereader
os.environ['PATH'] += ':/data/pacbio_assembly/AwesomeAssembler/DALIGNER'
#print os.popen("export").read()

Qvd = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
    'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D',
    'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
    'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
    'Y']
Qvv = range(len(Qvd))[::-1]

QVdict = dict(zip(Qvd,Qvv))


dbname = sys.argv[1]
lasname = sys.argv[2]
n = int(sys.argv[3])
path =  os.getcwd()+'/'
coveragename = path + 'coverage.txt'

aln = []

coveragefile = linereader.copen(coveragename)

coverage = coveragefile.getline(n)
cov = coverage.split()[2:]
covx = []
covy = []
for item in cov:
    data = item.split(',')
    covx.append(int(data[0]))
    covy.append(int(data[1]))

qv = list(util.get_QV(path+dbname, [n]))[0]
qx = []
qy = []
ts = int(sys.argv[5])
for i in range(len(qv)):
    qx.append(i*ts)
    qy.append(QVdict[qv[i]])

for item in util.get_alignments2(path+dbname,path+lasname,[n]):
    aln.append(item)
    
if (len(aln) == 0):
    sys.exit()
#print aln[0:5]

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

#print [len(item) for item in alns]
#print [item[0:3] for item in aln]

alns.sort(key = lambda x:min([item[3] for item in x]))

#size_chunk = num/grid_size
#for i in range(grid_size):
#    aln[i*size_chunk:min((i+1)*size_chunk, num)] = sorted(aln[i*size_chunk:min((i+1)*size_chunk, num)],key = lambda x: x[4]-x[3] ,reverse=True)

fig = plt.figure(figsize = (15,10))
plt.axes()
ax1 = plt.subplot2grid((6,6), (0, 0), colspan=6, rowspan=4)
ax2 = plt.subplot2grid((6,6), (4, 0), colspan=6, rowspan=1, sharex = ax1)
ax3 = plt.subplot2grid((6,6), (5, 0), colspan=6, rowspan=1, sharex = ax1)

#plt.gca().axes.get_yaxis().set_visible(False)
l = aln[0][5]
tip = l/200
ed = l/50
grid_size = 1.0
ax1.set_xlim(-2000,l+2000)
ax1.set_ylim(-5,num*grid_size)

points = [[0,0], [l,0], [l+tip,grid_size/4], [l,grid_size/2], [0,grid_size/2]]
#rectangle = plt.Rectangle((0, 0), l, 5, fc='r',ec = 'none')
polygon = plt.Polygon(points,fc = 'r', ec = 'none', alpha = 0.6)
ax1.add_patch(polygon)

dotted_line = plt.Line2D((0, 0), (0, num*grid_size ),ls='-.')               
ax1.add_line(dotted_line)

dotted_line2 = plt.Line2D((l, l), (0, num*grid_size ),ls='-.')               
ax1.add_line(dotted_line2)

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
        polygon.set_url('aln_svg' + str(item[2])+'.svg')
        ax1.add_patch(polygon)

        if points_end != []:
            polygon2 = plt.Polygon(points_end,fc = 'g', ec = 'none', alpha = 0.6)
            ax1.add_patch(polygon2)

        if points_start != []:
            polygon2 = plt.Polygon(points_start,fc = 'g', ec = 'none', alpha = 0.6)
            ax1.add_patch(polygon2)


ax2.plot(covx, covy)
ax3.plot(qx, qy)

plt.xlabel('position')
ax1.set_ylabel('pile-o-gram')
ax2.set_ylabel('coverage')
ax3.set_ylabel('i-qv')


plt.savefig(path + sys.argv[4] + '/aln_svg' + str(n)+ '.svg')