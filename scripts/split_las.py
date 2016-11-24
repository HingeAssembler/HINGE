#!/usr/bin/env python

import os
import argparse

ap = argparse.ArgumentParser(description="run LAsplit by splitting las into sizes of less than specified length")
ap.add_argument("las", help="path to las file to be split. assumed to be sorted.")
ap.add_argument("max_size", help="max size of any split file.", type=int, default=4, nargs='?')

args = ap.parse_args()

laspath = args.las
max_las_size = args.max_size


x = os.path.getsize(laspath)
num_divisions = (x/10**9)/max_las_size + 10
out_las_name = laspath.split('.las')[0]+'.# '

LAsplit_cmd = 'LAsplit -v '+out_las_name+ str(num_divisions) +' < ' + laspath
os.system(LAsplit_cmd)