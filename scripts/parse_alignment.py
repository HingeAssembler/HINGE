#!/usr/bin/python

import sys
import re


def parse_alignment(stream = sys.stdin):    
    with stream as f:
        for l in f:
            sub = re.sub('[\[\].x:<difstraep()]',' ',l.strip())    
            sub = re.sub(',','',sub)
            lst = sub.split()[:-1]
            if len(lst) == 9:
                yield [lst[2]] + map(int, lst[0:2] + lst[3:])