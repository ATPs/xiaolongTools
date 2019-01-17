import glob
import os
files = glob.glob('./potential_repeats/scaffold*')
print('number of files to combine', len(files))
results = [open(f).read() for f in files]
open('repeats.fa','w').write(''.join(results))