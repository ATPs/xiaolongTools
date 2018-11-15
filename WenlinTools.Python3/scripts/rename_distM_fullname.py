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
#import pandas as pd
#import datetime
from fullname_lib import get_names

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def xlsx2df(fn, sheet='Sheet1'):
    xl = pd.ExcelFile(fn)
    df = xl.parse(sheet)
    return df

def df2spdict(df, columns):
    adict = {}
    for index, row in df.iterrows():
        sp = row[columns[0]]
        sp = str(sp).split('-')[-1]
        line = [sp]
        for column in columns[1:]:
            each = row[column]
            if pd.isnull(each):
                continue
            else:
                if type(each) is datetime.datetime:
                    each = each.strftime('%Y-%m-%d')

                each = each.encode('ascii', 'ignore')
                line.append(each)

        adict[sp] = '\t'.join(line)
    return adict



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

def parse_spName(name, label):
    items = name.split('|')
    sp = items[0].split('_')[-1]
    ID = items[1]
    newName = []
    if 'coenia' in sp:
        newName.append('Coe')
    if 'grisea' in sp:
        newName.append('gri')
    if 'MEXICANspecies' in sp:
        newName.append('Msp')
    if 'neildi' in sp:
        newName.append('nld')
    if 'nigrosuffusa' in sp:
        newName.append('nig')
    if 'zonalis' in sp:
        newName.append('zol')

    if 'TX' in sp:
        newName.append('T')
    if 'hybrid' in sp:
        newName.append('H')

    if len(newName) == 0:
        newName.append(sp[:4])
    newName.insert(0, label)
    newName.append(ID)

    final = ''.join(newName)
    #final = final[:12]
    return final



if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    nameDict = get_names()


    info = []
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            info.append(line)
            continue
        items = line.strip().split()
        if len(items) == 1:
            info.append(line + '\n')
            continue
        sp = items[0].split('_')[0]
        if 'cp1' in line:
            label = 'A'
        else:
            label = 'B'
        name = parse_spName(nameDict[sp], label)
        print('nameCheck', name)
        try:
            line = '%s\t%s\n' % (name, '\t'.join(items[1:]))
        except:
            line = line + '\n'
        info.append(line)

    info = ''.join(info)

    dn = fn + '.renamed'
    cmn.write_file(info, dn)



