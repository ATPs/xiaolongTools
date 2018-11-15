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
        fn=sys.argv[1]
    except:
        print("Usage: *.py alist", file=sys.stderr)
        sys.exit()

    all_deflines = cmn.cmd2lines('grep ">" /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa|cut -d ">" -f 2')

    all_genus = set([each.split('_')[0] for each in all_deflines])

    for word in cmn.txt_read(fn).strip().split():
        if word in all_genus:
            print(word)


