import os
import random
import time
import sys

file_mzML = sys.argv[1]
file_index = sys.argv[2]

# file_mzML = '/gpfs/gpfs/project1/jx76-001/xc/MS/pride/mzML/PXD000004/40.mzML.gz'
# file_index = '/gpfs/gpfs/project1/jx76-001/xc/MS/proteins/20201201combined/20201202comet/20201201known.CONT.pep.join.idx'
# file_index = '/gpfs/gpfs/project1/jx76-001/xc/MS/proteins/20201201combined/20201202comet/20201201known.CONT.NEW.pep.join.idx'
tempfolder = '/home/scratch/'
tempfolder2 = '/gpfs/gpfs/project1/jx76-001/xc/MS/MSsearch/comet/temp/'
file_param = '/gpfs/gpfs/project1/jx76-001/xc/MS/proteins/20201201combined/comet.params.NoDecoySearch.NoCut.HCD'

basename_mzML = os.path.basename(file_mzML)
IDshort = basename_mzML.split('.')[0]
project_id = os.path.basename(os.path.dirname(file_mzML))
basename_index = os.path.basename(file_index)
outfolder_db = basename_index.rsplit('.',maxsplit=1)[0]
outfolder_pin = os.path.join('/gpfs/gpfs/project1/jx76-001/xc/MS/MSsearch/comet/', outfolder_db, project_id)
if not os.path.exists(outfolder_pin):
    time.sleep(random.random()*2)
    if not os.path.exists(outfolder_pin):
        os.makedirs(outfolder_pin)
outfile_pin = os.path.join(outfolder_pin, IDshort + '.pin.gz')
if os.path.exists(outfile_pin) and os.path.getsize(outfile_pin) > 10000:
    print('seems finished', IDshort)
    exit()

# move index file to tempfolder if it does not exist
file_index2 = os.path.join(tempfolder, basename_index)
time.sleep(random.random()*2)
if not os.path.exists(file_index2):
    print(f'cp index file to {tempfolder}')
    os.system(f'cp {file_index} {tempfolder}')
else:
    for i in range(300):
        time.sleep(1)
        if os.path.getsize(file_index) == os.path.getsize(file_index2):
            print(f'index file already in {tempfolder}')
            break
    else:
        print('too slow copy the index file, please check')
        exit()

# cp mzML file to tempfolder. if mzML file greater than 5GB, decompress it
mzML_size = os.path.getsize(file_mzML)
if mzML_size > 5e9:
    tempWork = os.path.join(tempfolder2, basename_index.rsplit('.',maxsplit=1)[0] + '_' + IDshort)# working folder
else:
    tempWork = os.path.join(tempfolder, basename_index.rsplit('.',maxsplit=1)[0] + '_' + IDshort)# working folder
if not os.path.exists(tempWork):
    os.makedirs(tempWork)

if mzML_size > 5e9:
    os.system(f'ln -s {file_mzML} {tempWork}/')
else:
    os.system(f'cp {file_mzML} {tempWork}/')

# run comet
cmd = f'cd {tempWork} && comet -P{file_param} -D{file_index2} {basename_mzML} && pigz -c {IDshort}.pin > {outfile_pin}; rm -rf {tempWork}'
os.system(cmd)