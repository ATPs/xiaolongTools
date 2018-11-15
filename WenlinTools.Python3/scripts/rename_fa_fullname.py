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
#import ete3
from fullname_lib import get_names
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def old_get_names():
    adict = {}
    #1. firstly get data download manually
    fn1 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/All-plates.xlsx'
    df = xlsx2df(fn1)
    taken_columns = [
            'DNA number',
            'Taxon name',
            'Type status',
            'Sex',
            'Country',
            'State/Province',
            'County/Region',
            'Date',
            ]
    df = df.loc[:, taken_columns]

    adict.update(df2spdict(df, taken_columns))

    fn2 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/BugsInRNAlater-2015-04-11.xlsx'
    df = xlsx2df(fn2)
    taken_columns = [
            'Number',
            'Species name',
            'Sex',
            'State/Country',
            'County',
            'Date'
            ]
    df = df.loc[:, taken_columns]

    adict.update(df2spdict(df, taken_columns))

    return adict


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

    nameDict = get_names()

    #t = ete3.Tree(cmn.txt_read(fn).replace('[&U]', ''))

    appear = {}
    table = []
    lines = cmn.file2lines(fn)
    info = []
    #check appearance
    for line in lines:
        sp = line[1:].split('_')[0].replace('flt', '').replace('.Lerema.out', '')
        if True:
            if sp in appear:
                appear[sp] += 1
            else:
                appear[sp] = 1

    multiple_sps = [sp for sp in appear
            if appear[sp] > 1]

    appear = {}
    for line in lines:
        if line[0] =='>':
            sp = line[1:].split('_')[0].replace('flt', '').replace('.fa', '').replace('.Lerema.out', '')
            try:
                #newline = line.replace(sp, nameDict[sp].replace('\t', '_').replace('  ', ' ').replace(' ', '_'))
                newline = '>%s' % (nameDict[sp].replace('\t', '_').replace('  ', ' ').replace(' ', '_'))
            except:
                print('can not find name for %s' % sp)
                newline = line

            if sp in appear:
                appear[sp] += 1
                #newline += '_cp%s' % appear[sp]
            else:
                appear[sp] = 1
                #if sp in multiple_sps:
                #    newline += '_cp%s' % appear[sp]
            info.append(newline)
        else:
            info.append(line)

    info.append('')
    info = '\n'.join(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)
    cmn.write_file(''.join(table), dn + '.nameTable')



