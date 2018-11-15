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
        fqlist, fout = sys.argv[1:]
    except:
        print("Usage: *.py fqlist fout", file=sys.stderr)
        sys.exit()


    badones = set([])

    with open(fout) as fp:
        for line in fp:
            items = line.strip().split()
            sseqid = items[1]
            pident = float(items[2])
            evalue = float(items[3])
            if pident == 100:
                badones.add(sseqid)

            if evalue <= 0.01:
                badones.add(sseqid)



    for fq in cmn.file2lines(fqlist):
        dn = 'noWolb_%s' % cmn.lastName(fq)
        with open(fq) as fp, open(dn, 'w') as dp:
            for i, line in enumerate(fp):
                if i % 4 == 0:
                    ID = line.strip()[1:].replace('/', '_').replace(' ', '_')
                    if ID in badones:
                        isGood = False
                    else:
                        isGood = True

                if isGood:
                    dp.write(line)
