import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

IDs = set([])

fn = sys.argv[1]

with open(fn) as fp:
    for i, line in enumerate(fp):
        if i % 4 == 0:
            ID = ':'.join(line.split(':')[:3])
            IDs.add(ID)


print('\n'.join(IDs))

        
