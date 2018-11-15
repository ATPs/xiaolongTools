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
        fns=sys.argv[1:]
    except:
        print("Usage: *.py read_filelist", file=sys.stderr)
        sys.exit()

    #read in data
    fdict = 'blastBySp.dict.pkl'
    sp_dict = cmn.pickle_read(fdict)

    #get the read ID, and the exon of it
    good_IDs = {}
    for sp in sp_dict:
        lines = sp_dict[sp]
        for line in lines:
            readID = line.split()[2]
            print('readID', readID)
            good_IDs[readID] = sp


    #get the reads and split them into exons
    #fns = cmn.getid(fn)

    rdict = {}
    for fn in fns:
        print('parsing ' + fn)
        with open(fn) as fp:
            for i, line in enumerate(fp):
                if i % 4 == 0:
                    #record = []
                    ID = line.strip().replace(' ', '_')
                    #print 'checkID', ID
                    try:
                        sp = good_IDs[ID]
                        isGood = True
                    except KeyError:
                        isGood = False

                #record.append(line)

                if i % 4 == 1:
                    #record = ''.join(record)
                    seq = line.strip()
                    if isGood:
                        if sp not in rdict:
                            rdict[sp] = {}
                        rdict[sp][ID] = seq

    dn = 'readsBySp.dict.pkl'
    cmn.pickle_write(rdict, dn)



