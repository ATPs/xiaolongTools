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
import barcode_processing as bp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py fa", file=sys.stderr)
        sys.exit()


    seqDict = read_fa(fn)

    conn = bp.lock_database()
    c = conn.cursor()
    hasError = False
    for name in seqDict:
        genus, sp = bp.find_genus_and_species(name)
        if '|' not in name and (name.split('_') <= 2):
            hasError = True
            print('Error: please put identifier to the species name. for example, "Junonia_coenia" is not acceptable but "Junonia_coenia_3935" is accepted')
            print('identifier could be anything, such as date or location')
        #sampleID, CSid = bp.find_IDs(name)
        #if sampleID != 'NA' or CSid != 'NA':
        #    print 'found %s %s in %s' % (sampleID, CSid, name)
        c.execute('select * from barcode_fasta where full_name = \'%s\'' % name)
        lines = c.fetchall()
        if len(lines) != 0:
            hasError = True
            print('Error! the name "%s" exists in the current dataset! please change it to a new name or add some affix' % name)
        if genus == 'NA':
            hasError = True
            print('can not find genus for (%s, %s) %s' % (genus, sp, name))
        seq = seqDict[name]
        if len(seq) != 658:
            hasError = True
            print('Error! this sequence length is not 658: %s %s' % (name, len(seq)))

    conn.close()
    if hasError:
        print('\nPlease fix those errors before proceeding')
    else:
        print('\nso far so good!')
