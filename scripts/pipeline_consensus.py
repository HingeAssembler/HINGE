#!/usr/bin/python

import sys
import os
import subprocess


if len(sys.argv) >= 2:
	bact_id = sys.argv[1]


ini_path = 'nominal.ini'
if len(sys.argv) >= 3:
    ini_path = sys.argv[2]


run_identifier = 'A'
if len(sys.argv) >= 4:
    run_identifier = sys.argv[3]

graphml_file = bact_id+run_identifier+'.G2.graphml'


# This is used to start the pipeline in the middle
st_point = 0
if len(sys.argv) >= 5:
	st_point = int(sys.argv[4])

base_path = './'

if st_point <= 1:
    draft_path_cmd = 'get_draft_path.py '+base_path+' '+ bact_id+' '+graphml_file
    print '1: '+draft_path_cmd
    subprocess.check_output(draft_path_cmd,cwd=base_path, shell=True)


if st_point <= 2:
    draft_assembly_cmd = 'draft_assembly --db '+bact_id+' --las '+bact_id+'.las --prefix '+bact_id+' --config '+ini_path+' --out '+bact_id+'.draft'
    print '2: '+draft_assembly_cmd
    subprocess.check_output(draft_assembly_cmd,cwd=base_path, shell=True)
  

if st_point <= 3:
    corr_head_cmd = 'correct_head.py '+bact_id+'.draft.fasta '+bact_id+'.draft.pb.fasta draft_map.txt'
    print '3: '+corr_head_cmd
    subprocess.check_output(corr_head_cmd,cwd=base_path, shell=True)


if st_point <= 4:
    subprocess.call("rm -f draft.db",shell=True,cwd=base_path)
    fasta2DB_cmd = "fasta2DB draft "+base_path+bact_id+'.draft.pb.fasta'
    print '4: '+fasta2DB_cmd
    subprocess.check_output(fasta2DB_cmd.split(),cwd=base_path)

if st_point <= 5:
    subprocess.call("rm -f draft.*.las",shell=True,cwd=base_path)
    mapper_cmd = "HPCmapper draft "+bact_id
    print '5: '+mapper_cmd
    subprocess.call(mapper_cmd.split(),stdout=open(base_path+'draft_consensus.sh','w') , cwd=base_path)

if st_point <= 6:
    mapper_shell_cmd = "csh -v draft_consensus.sh"
    print '6: '+mapper_shell_cmd
    subprocess.check_output(mapper_shell_cmd.split(), cwd=base_path)

if st_point <=7:
    # remove_cmd = 'rm -f nonrevcompdraft.'+bact_id+'.*.las'
    # subprocess.call(remove_cmd,shell=True,cwd=base_path)
    LAmerge_cmd = "LAmerge draft."+bact_id+".las "+'draft.'+bact_id+'.[0-9].las'
    print '7: '+LAmerge_cmd
    subprocess.check_output(LAmerge_cmd,cwd=base_path,shell=True)

if st_point <= 8:
    consensus_cmd = 'consensus draft '+bact_id+' draft.'+bact_id+'.las '+bact_id+'.consensus.fasta '+ini_path
    print '8: '+consensus_cmd
    subprocess.check_output(consensus_cmd,cwd=base_path,shell=True)
    

if st_point <= 9:
    gfa_cmd =  'get_consensus_gfa.py '+base_path+ ' '+ bact_id+ ' '+graphml_file+' '+bact_id+'.consensus.fasta' 
    print '9: '+gfa_cmd
    subprocess.check_output(gfa_cmd,cwd=base_path,shell=True)




