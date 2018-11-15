import pandas as pd
import datetime

tableDir = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/downloaded_excels'
#tableDir = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables'

#old path: /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/downloaded_excels
#ming path: /work/biophysics/mtang/SNP_calling/scripts/data/nameTables
def xlsx2df(fn,sheet='Sheet1'):
    xl = pd.ExcelFile(fn)
    df = xl.parse(sheet)
    return df

def df2spdict(df, columns):
    adict = {}
    Sdict = {}
    IDlabels = ['SRNP', 'CSUPOBK', 'CSU-CPG-LEP']
    for index, row in df.iterrows():
        sp = row[columns[0]]

        try:
            sp = str(sp).strip().replace('LEP-', 'LEP').replace('NVG-', '').replace('11-BOA-', '').replace('YPM-ENT-', '').strip()
        except:
            continue
        line = [sp]
        for eachCell in row:
            try:
                eachCell = str(eachCell)
            except:
                continue

            if any([IDlabel in eachCell
                for IDlabel in IDlabels]):
                Sdict[sp] = eachCell.replace('?', '')


            #CSUPOBK or CSU-CPG-LEP
                Sdict[sp] = eachCell.replace('?', '')

        for column in columns[1:]:
            each = row[column]
            if pd.isnull(each):
                continue
            else:
                if type(each) is datetime.datetime:
                    each = each.strftime('%Y-%m-%d')

                try:
                    each = each.encode('ascii', 'ignore')
                except:
                    continue
                each = each.replace('"', '')
                line.append(each)

        adict[sp] = '\t'.join(line)

    #print Sdict
    return Sdict, adict



def get_SRNPnumber():
    adict = {}
    SRNP_dict = {}
    #1. firstly get data download manually
    #fn1 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/All-plates.xlsx'
    fn1 = '%s/All-plates.xlsx' % tableDir
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
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn2 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/BugsInRNAlater-2015-04-11.xlsx'
    fn2 = '%s/BugsInRNAlater-2015-04-11.xlsx' % tableDir
    df = xlsx2df(fn2)
    taken_columns = [
            'DNA Number',
            'Taxon name',
            'Sex',
            'State',
            'County/Region',
            'Date'
            ]
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn3 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/All-sequenced.xlsx'
    fn3 = '%s/All-sequenced.xlsx' % tableDir
    print('all-sequenced')
    df = xlsx2df(fn3, 'All')
    taken_columns = [
            'DNA number',
            'Taxon name',
            'Sex',
            'Country',
            'State/Province',
            'County/Region',
            'Date',
            ]
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn4 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/Junonia-all-sequenced.xlsx'
    fn4 = '%s/Junonia-all-sequenced.xlsx' % tableDir
    df = xlsx2df(fn4)
    taken_columns = [
            'Number',
            'Treename'
            ]
    df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    return SRNP_dict


def get_names():
    adict = {}
    SRNP_dict = {}
    #1. firstly get data download manually
    #fn3 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/All-sequenced.xlsx'
    #print 'all-sequenced'
    #df = xlsx2df(fn3, 'All')
    #taken_columns = [
    #        'DNA number',
    #        'Taxon name',
    #        'Sex',
    #        'Country',
    #        'State/Province',
    #        'County/Region',
    #        'Date',
    #        ]
    #df = df.loc[:, taken_columns]

    #S_dict, subdict = df2spdict(df, taken_columns)
    #adict.update(subdict)
    #SRNP_dict.update(S_dict)

    #fn1 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/All-plates.xlsx'
    fn1 = '%s/All-plates.xlsx' % tableDir
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
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn2 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/BugsInRNAlater-2015-04-11.xlsx'
    fn2 = '%s/BugsInRNAlater-2015-04-11.xlsx' % tableDir
    df = xlsx2df(fn2)
    taken_columns = [
            'DNA Number',
            'Taxon name',
            'Type status',
            'Sex',
            'Country',
            'State',
            'County/Region',
            'Date'
            ]
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn4 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/Junonia-all-sequenced.xlsx'
    fn4 = '%s/Junonia-all-sequenced.xlsx' % tableDir
    df = xlsx2df(fn4)
    taken_columns = [
            'Number',
            'Treename'
            ]
    df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #add SRNP number into adict
    for sample in SRNP_dict:
        adict[sample] += '\t' + SRNP_dict[sample]

    return adict


def get_names_4barcode():
    adict = {}
    SRNP_dict = {}
    #1. firstly get data download manually
    #fn1 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/All-plates.xlsx'
    fn1 = '%s/All-plates.xlsx' % tableDir
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
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn2 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/BugsInRNAlater-2015-04-11.xlsx'
    fn2 = '%s/BugsInRNAlater-2015-04-11.xlsx' % tableDir
    df = xlsx2df(fn2)
    taken_columns = [
            'DNA Number',
            'Taxon name',
            'Sex',
            'State',
            'County/Region',
            'Date'
            ]
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    #fn3 = '/work/00412/mtang/sequencing/scripts/downloaded_excels/All-sequenced.xlsx'
    #print 'all-sequenced'
    #df = xlsx2df(fn3, 'All')
    #taken_columns = [
            #'DNA number',
            #'Taxon name',
            #'Sex',
            #'Country',
            #'State/Province',
            #'County/Region',
            #'Date',
            #]
    #df = df.loc[:, taken_columns]

    #S_dict, subdict = df2spdict(df, taken_columns)
    #adict.update(subdict)
    #SRNP_dict.update(S_dict)

    #fn4 = '/work/00412/mtang/sequencing/scripts/downloaded_excels/Junonia-all-sequenced.xlsx'
    #df = xlsx2df(fn4)
    #taken_columns = [
            #'Number',
            #'Treename'
            #]
    #df = df.loc[:, taken_columns]

    #S_dict, subdict = df2spdict(df, taken_columns)
    #adict.update(subdict)
    #SRNP_dict.update(S_dict)

    #add SRNP number into adict
    for sample in SRNP_dict:
        adict[sample] += '\t' + SRNP_dict[sample]

    return adict


def get_IDlist():
    #fn3 = '/work/biophysics/mtang/SNP_calling/scripts/data/nameTables/All-sequenced.xlsx'
    fn3 = '%s/All-sequenced.xlsx' % tableDir
    #print 'all-sequenced'
    df = xlsx2df(fn3, 'All')
    taken_columns = [
            'DNA number',
            ]
    #df = df.loc[:, taken_columns]
    thelist = df[taken_columns[0]].tolist()
    new = []
    for ID in thelist:
        try:
            ID = str(ID)
        except:
            continue
        new.append(ID)

    return new

    #return list(df)
    #return df.tolist()
