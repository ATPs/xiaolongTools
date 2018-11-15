import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

template = cmn.txt_read('/work/biophysics/mtang/SNP_calling/scripts/templates/filter_gap.template')

input = sys.argv[1]

info = template.replace('[INPUT]', input)

dn = 'filter.job'
cmn.write_file(info, dn)

print('please submit %s to the queue' % dn)

