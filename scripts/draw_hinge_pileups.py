#!/usr/bin/python

import sys
import os





def draw_hinges(hinge_file):

    with open (hinge_file) as f:

    	bact_name = hinge_file.split('.')[0]

        for lines in f:
            lines1=lines.split()
            node_id = int(lines1[0]) + 1

            os.system("python ~/AwesomeAssembler/scripts/draw2.py "+bact_name+" "+bact_name+".las "+str(node_id)+" sample 100")
  

if __name__ == "__main__":   

    draw_hinges(sys.argv[1])





