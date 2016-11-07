#!/usr/bin/env python

import sys
import os



def reverse_complement(bases):
    rev_comp={'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return ''.join(map(lambda x :rev_comp[x],bases[::-1]))

def run(multifasta_path,intermediate_repeat_file_path, gt_file_path, gt_annotated_file_path):

	##Read chromosomes from multifasta file
	chrom={}
	cur_chrom=''
	start=True
	chr_num=0
	with open(multifasta_path,'r') as f:
	    for lines in f:
	        if lines[0]=='>':
	            
	            if start:
	                start=False
	                chr_num=int(lines.split()[0][1:])-1
	                print 'detect chr '+ str(chr_num)
	            else:
	                chrom[chr_num]=cur_chrom
	                print len(cur_chrom)
	                chr_num=int(lines.split()[0][1:])-1
	                print 'detect chr '+ str(chr_num)
	            cur_chrom=''
	        else:
	            cur_chrom+=lines.strip()
	print len(cur_chrom)            
	chrom[chr_num]=cur_chrom

	##Run mummer to get repeats
	mummer_cmd='mummer -maxmatch -b -c -l 1000 -L '+multifasta_path+' '+multifasta_path +' > '+intermediate_repeat_file_path 
	os.system(mummer_cmd)


	#Put repeats discovered by mummer in right form
	chr_num=0
	chr_repeats={}
	rev_com=False
	with open (intermediate_repeat_file_path) as f:
	    for line in f:
	        if line[0]=='>':
	            line1=line.strip().split()
	            chr_num=int(line1[1])-1
	            #print len(line1)
	            if len(line1)==6:
	                rev_com=True
	                chr_len=int(line1[5])
	            else:
	                rev_com=False
	                chr_len=int(line1[4])
	        else:
	            line1=line.strip().split()
	            chr2_num=int(line1[0])-1
	            chr2_start=int(line1[1])-1
	            chr1_start=int(line1[2])-1
	            rep_len=int(line1[3])
	            if not rev_com:
	                if not (chrom[chr_num][chr1_start:chr1_start+rep_len]
	                        ==chrom[chr2_num][chr2_start:chr2_start+rep_len]):
	                    print chr_num+1,line
	                if chr1_start==0 and rep_len==chr_len:
	                    continue
	                chr_repeats.setdefault(chr_num,[]).append((chr1_start,chr1_start+rep_len))
	            else:
	                if not (chrom[chr_num][chr1_start-rep_len+1:chr1_start+1]
	                        == reverse_complement(chrom[chr2_num][chr2_start:chr2_start+rep_len])):
	                    print chr_num+1,line,rev_com
	                chr_repeats.setdefault(chr_num,[]).append((chr1_start-rep_len+1,chr1_start+1))

	#Go through gt file and annotate reads that intersect with repeats.
	with open(gt_file_path) as f:
	    with open(gt_annotated_file_path,'w') as g:
	        for line in f:
	            line1=line.split()
	            cr=int(line1[1])
	            rd_st=int(line1[2])
	            rd_end=int(line1[3])
	            is_repeat=0
	            for tup in chr_repeats[cr]:
	                if ((rd_st >= tup[0] and rd_st <= tup[1]) or
	                    (rd_end >= tup[0] and rd_end <= tup[1])):
	                    is_repeat=1
	            line2=line.strip()+"\t"+str(is_repeat)+"\n"
	            g.write(line2)

if __name__ == '__main__':
	multifasta_path=sys.argv[1]
	gt_file_path=sys.argv[2]
	gt_annotated_file_path=sys.argv[3]
	intermediate_repeat_file_path='./repeats_discovered.txt'
	if len(sys.argv) > 4:
		intermediate_repeat_file_path=sys.argv[4]
	run(multifasta_path,intermediate_repeat_file_path, gt_file_path, gt_annotated_file_path)ß
