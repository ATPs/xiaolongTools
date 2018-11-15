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
import ete3
import pandas as pd
import datetime

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def xlsx2df(fn,sheet='Sheet1'):
    xl = pd.ExcelFile(fn)
    df = xl.parse(sheet)
    return df

def df2spdict(df, columns):
    adict = {}
    for index, row in df.iterrows():
        sp = row[columns[0]]
        sp = str(sp).split('-')[-1]
        line = []
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



def get_names():
    adict = {}
    #1. firstly get data download manually
    fn4 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/Junonia-all-sequenced.xlsx'
    df = xlsx2df(fn4)
    taken_columns = [
            'Number',
            'Treename'
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

    t = ete3.Tree(cmn.txt_read(fn).replace('[&U]', ''))

    appear = {}
    table = []
    for node in t:
        name = node.name
        sp = name.split('_')[0].split('.')[0]
        if sp not in appear:
            appear[sp] = 1
            table.append('%s\n' % ( nameDict[sp]))
        else:
            appear[sp] += 1

        new_name = '%s_%s' % (nameDict[sp], appear[sp])
        #new_name = name.replace(sp, nameDict[sp])
        node.name = new_name

    info = t.write()
    print(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)
    cmn.write_file(''.join(table), dn + '.nameTable')



