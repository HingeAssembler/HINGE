import json
import os
import sys
import subprocess

base_dir = '/data/pacbio_assembly/pb_data/NCTC/'
bact_dict = json.load(open(base_dir+'NCTC.json'))

#bacterium_of_interest='NCTC7972'

bacterium_of_interest=sys.argv[1]

if len(sys.argv) > 2:
	bact_dict=sys.argv[2]

bact_name="_".join(bact_dict[bacterium_of_interest]['Species'][0].split())

cmd_base = 'ascp -QT -l 1000m -i /data/pacbio_assembly/pb_data/asperaweb_id_dsa.openssh era-fasp@fasp.ega.ebi.ac.uk:vol1/'
dest_dir = base_dir+bacterium_of_interest+'/'

os.system('mkdir -p '+dest_dir)

for run, file_list in bact_dict[bacterium_of_interest]['file_paths'].items():
    for file_path in  file_list:
        cmd = cmd_base+file_path+' '+dest_dir
        print cmd
        os.system(cmd)

dest_fasta_name = dest_dir+bact_name

dextract_cmd = 'dextract -o'+dest_fasta_name

bax_files = [x for x in os.listdir(dest_dir) if x.endswith('.bax.h5')]

for bax_file in bax_files:
	dextract_cmd +=  " " + dest_dir+bax_file

print dextract_cmd
dextract_cmd += 'lkasjdfla.bax.h5'
subprocess.call(dextract_cmd.split())

print 'done'



