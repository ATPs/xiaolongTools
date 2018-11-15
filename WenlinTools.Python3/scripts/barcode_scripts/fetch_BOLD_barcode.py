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
        print("Usage: *.py link.file", file=sys.stderr)
        sys.exit()


    dn = 'retrieved_barcodes.fa'

    dp = open(dn, 'w')
    for link in cmn.file2lines(fn):
        if link[0] == '#':
            continue
        print('processing ' + link)
        info = cmn.link2info(link)
        seq = cmn.find_between(info, "generateBarcode ('#barcodeImg_", "');").split('\'')[-1]
        takeName = False
        takeSp = False
        for line in info.split('\n'):
            if 'Sequence ID' in line:
                takeName = True
                continue
            if '<td>Species:</td>' in line:
                takeSp = True
                continue

            if takeName:
                takeName = False
                #<td style="width:160px;">ANICE505-10.COI-5P</td>
                name = cmn.find_between(line, '>', '<').split('.COI')[0]
            if takeSp:
                takeSp = False
                sp = cmn.find_between(line, '<em>', '</em>').replace(' ', '_')
        defline = '%s_%s' % (sp, name)
        print(defline)
        fasta = '>%s\n%s\n' % (defline, seq)
        dp.write(fasta)

    dp.close()


