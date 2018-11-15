#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    f_table = '/project/biophysics/Nick_lab/wli/sequencing/scripts/name_table'

    nameDict = {}
    for line in cmn.file2lines(f_table):
        items = line.strip().split()
        if len(items) == 0:
            continue
        label = items[0]
        name = '_'.join(items[1:])
        nameDict[label] = name.replace('-', '_')

    print(list(nameDict.keys()))

    info = cmn.txt_read(fn)
    hasDash = ('_' in info)
    for label in nameDict:
        name = nameDict[label]
        #data_label = '%s_' % label
        #if label.isdigit():
        if True:
            if hasDash:
                data_label = '%s_' % label
                new_name = '%s_%s_' % (label, name)
            else:
                data_label = '%s:' % label
                new_name = '%s_%s:' % (label, name)

        else:
            data_label = '%s' % label
            #new_name = '%s_%s' % (label, name)
            new_name = '%s_%s' % (label, name)
        info = info.replace(data_label, new_name)

    info = info.replace('_snp.codeVcf', '')
    info = info.replace('.codeVcf', '')
    info = info.replace('_father', '_copy1')
    info = info.replace('_mother', '_copy2')
    print(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)



