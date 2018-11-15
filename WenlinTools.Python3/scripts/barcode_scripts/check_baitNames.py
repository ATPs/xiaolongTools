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

    goodIDlist = set([])
    allIDdict = {}
    goodLines = set([])
    for line in cmn.file2lines(fn):
        items = line.strip().split()
        if len(items) != 2:
            print('Error! please only put sampleID and genus name in one line!')
            hasError = True
            continue
        else:
            sp, genus = items
            try:
                allIDdict[sp].append(genus)
            except KeyError:
                allIDdict[sp] = [genus]

            if genus in all_genus:
                #print 'Error! genus %s didn\'t have barcode recorded!' % genus
                #hasError = True
                goodIDlist.add(sp)
                goodLines.add(line)


    hasError = False
    #if not hasError:
    #    print 'so far so good!'
    for sp in allIDdict:
        if sp not in goodIDlist:
            hasError = True
            print('Error: all the genues for this sample (%s) is not present:' % sp)
            print('Error: genus list: %s' % ','.join(list(allIDdict.values())))

    if not hasError:
        print('so far so good')

    #modify the file to remove bad ones
    fbad = 'removed_badOnes.list'
    ftmp = cmn.lastName(fn) + '.tmp'
    cmn.run('cp %s %s' % (fn, ftmp))

    with open(ftmp) as fp, open(fbad, 'a') as dpBad, open(fn, 'w') as dp:
        for line in fp:
            if line.strip() in goodLines:
                dp.write(line)
            else:
                dpBad.write(line)

