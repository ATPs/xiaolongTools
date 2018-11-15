#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os
#import pandas as pd
#import datetime
from fullname_lib import get_names_4barcode

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0].split()[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


def read_autoTable(fn):
    seqDict = read_fa(fn)

    new = {}
    for name in seqDict:
        seq = seqDict[name]
        name = '_'.join(name.split('_')[1:])
        new[name] = seq
    return new


def attempt_to_find_genus_by_abundence(ID, fqlist):
    tmpdir = 'tmp_%s' % ID
    cmn.mkdir(tmpdir)
    os.chdir(tmpdir)

    cmn.write_lines(fqlist, 'fqlist')
    cmd = '/work/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/auto_rebait.py fqlist'
    cmn.run(cmd)

    dn = 'picked_bait.txt'
    if cmn.filexist(dn):
        genus = cmn.txt_read(dn).strip().split('_')[0].split()[0]
    else:
        genus = None
    os.chdir('..')
    cmn.run('cp %s/mapping_stat.info tmpStat/%s_mapping_stat.info' % (tmpdir, ID))
    cmn.run('rm -r %s ' % tmpdir)
    return genus



if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py fqlist", file=sys.stderr)
        sys.exit()

    cmn.mkdir('tmpStat')

    IDlist = set([])
    fq_groups = {}
    for line in cmn.file2lines(fn):
        Id = cmn.lastName(line).split('_')[0]
        Id = Id.replace('NVG-', '').replace('11-BOA-','').replace('LEP-', 'LEP')
        IDlist.add(Id)
        fq = os.path.abspath(line)
        try:
            fq_groups[Id].append(fq)
        except KeyError:
            fq_groups[Id] = [fq]

    nameDict = get_names_4barcode()

    fall = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #fall = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes.fasta'
    seqDict = read_fa(fall)
    fadd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/addedBaits_fromPipeline.fa'
    if cmn.filexist(fadd):
        seqDict.update(read_fa(fadd))
    ftable = '/archive/biophysics/Nick_lab/wli/archive/barcodes/auto_tables/verified_barcodes.fa'
    seqDict.update(read_autoTable(ftable))
    all_genus = set([name.split('_')[0].lower() for name in seqDict])

    info = []
    missing = []
    notFound = []

    for sp in IDlist:
        try:
            fullname = nameDict[sp]
            genus = fullname.split()[1]
            if genus.lower() not in all_genus:
                notFound.append(sp)

        except KeyError:
            missing.append(sp)
            continue

    if len(notFound) != 0:
        try:
            #if input such info, then don't need to run abundance check
            corrected_info = cmn.file2lines(sys.argv[2])
            resolvedIDs = [each.strip().split()[0] for each in corrected_info]
            notFound = set(notFound) - set(resolvedIDs)
        except:
            resolved = []
            resolvedIDs = []
            print('attempting to find the species by abundon...')
            for ID in notFound:
                fqlist = fq_groups[ID]
                genus = attempt_to_find_genus_by_abundence(ID, fqlist)
                if genus != None:
                    print('found %s for %s' % (genus, ID))
                    resolved.append('%s\t%s\n' % (ID, genus))
                    resolvedIDs.append(ID)

            if resolved != []:
                cmn.write_file(''.join(resolved), 'resolved_genus_by_abundance.txt')
            notFound = set(notFound) - set(resolvedIDs)


    print('\n\n\n\n')
    print('###############################################################')
    print('##################begin reports below##########################')
    print('###############################################################')
    if len(missing) == 0 and len(notFound) == 0:
        print('GoodNews: Every sample has a reference bait in their genus')

    if len(missing) != 0:
        print('\nthe following IDs are missing in all-sequenced table (update table and try again?):')
        print('\n'.join(missing))

    if len(notFound) != 0:
        print('\nPlease ask Nick for the bait of the following IDs:')
        for sp in notFound:
            print(nameDict[sp])
