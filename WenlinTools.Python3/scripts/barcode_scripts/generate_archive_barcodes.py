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
import re
import datetime
import pandas as pd
import xlsxwriter
import ete3
import barcode_processing as bp
gapN_dict = {}


def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict, len(seq)

def xlsx2df(fn,sheet='Sheet1'):
    print('reading table %s ...' % fn)
    xl = pd.ExcelFile(fn)
    df = xl.parse(sheet)
    return df

def df2spdict(df, columns):
    adict = {}
    Sdict = {}
    for index, row in df.iterrows():
        sp = row[columns[0]]
        sp = str(sp).replace('LEP-', 'LEP').split('-')[-1]
        line = [sp]
        for eachCell in row:
            try:
                eachCell = str(eachCell)
            except:
                continue

            if 'SRNP' in eachCell:
                Sdict[sp] = eachCell

        for column in columns[1:]:
            #TODO: how placeholder works here?
            try:
                each = row[column]
            except KeyError:
                line.append('')
                continue

            if pd.isnull(each):
                line.append('')
            else:
                if type(each) is datetime.datetime:
                    each = each.strftime('%Y-%m-%d')

                try:
                    each = str(each).encode('ascii', 'ignore')
                except:
                    tmp = []
                    for char in each:
                        try:
                            tmp += str(char)
                        except:
                            continue
                    each = ''.join(tmp)
                line.append(each)

        adict[sp] = line

    return Sdict, adict



def gather_spInfo():
    adict = {}
    SRNP_dict = {}
    header = 'sampleID Taxon_name Type_status Sex Country State County Date SRNP'.split()
    #1. firstly get data download manually
    #fn3 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/All-sequenced.xlsx'
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
    #df = df.loc[:, taken_columns]

    S_dict, subdict = df2spdict(df, taken_columns)
    adict.update(subdict)
    SRNP_dict.update(S_dict)

    fn2 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/BugsInRNAlater-2015-04-11.xlsx'
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

    #fn4 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/Junonia-all-sequenced.xlsx'
    #df = xlsx2df(fn4)
    #taken_columns = [
    #        'Number',
    #        'Treename'
    #        ]
    #df = df.loc[:, taken_columns]

    #S_dict, subdict = df2spdict(df, taken_columns)
    #adict.update(subdict)
    #SRNP_dict.update(S_dict)

    #add SRNP number into adict
    for sample in list(adict.keys()):
        try:
            adict[sample].append(SRNP_dict[sample])
        except:
            adict[sample].append('')

    return adict, header


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_fasta(lines):
    namelist = []
    seqDict = {}
    seqs = []
    name = None
    for line in lines:
        line = line.strip()
        if line[0] == '>':
            if seqs != []:
                seqDict[name] = ''.join(seqs)
                seqs = []

            name = line[1:]
            name = parse_name(name)
            namelist.append(name)

        else:
            seq = line
            seqs.append(seq)

    seqDict[name] = ''.join(seqs)

    for name in list(seqDict.keys()):
        seqDict[name] = ''.join(seqDict[name]).replace('-', 'N')

    return namelist, seqDict


def parse_name(name):
    name = name.replace('[new]', '').replace('[old]', '')
    return name


def parse_clw(lines):
    namelist = []
    seqDict = {}
    for line in lines:
        name, seq = line.strip().split()
        name = parse_name(name)
        namelist.append(name)
        seqDict[name] = seq.replace('-', 'N')
    return namelist, seqDict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def readin_barcode_file(lines):
    if any([line.strip()[0] == '>'
        for line in lines
        if line.strip() != '']):
        alist, seqDict = parse_fasta(lines)
    else:
        alist, seqDict = parse_clw(lines)

    return alist, seqDict


def tell_isPaired(names, fn):
    if cmn.lastName(fn) == 'unconfident_barcodes_paired.fa':
        return 'evenLine'

    #use the (assembled) label
    assembled_lines = [name for name in names
            if '(assembled)' in name]
    Nass = len(assembled_lines)
    print('checkLog', Nass, len(names))

    if 0 < Nass <= (len(names) / 2):
        return 'assembled'#is paired

    #check number for pairing
    Nass = 0
    for name in names:
        if name[0].isdigit():
            Nass += 1

    if Nass <= (len(names) / 2):
        return 'evenLine'

    if 'pair' in fn:
        return 'evenLine'

    return False


def tell_pair_sequence(pairType, i, name):
    #return True if this sequence is paired and NOT the orignal read
    #default would be False (original read)
    if pairType == False:
        return False
    else:
        if pairType == 'assembled':
            if '(assembled)' not in name:
                return True # paired read

        if pairType == 'evenLine':
            if i % 2 == 1:
                return True # paired read

        return False

def name2sp(name):
    labels = ['[new]', '[old]', '[unverified]', '[notSure]', '(assembled)']
    sp = name.strip().split('_')[0]
    for label in labels:
        sp = sp.replace(label, '')
    return sp

def remove_labels(name):
    labels = ['[new]', '[old]', '[unverified]', '[notSure]',]
    for label in labels:
        name = name.replace(label, '')
    return name


def print4sort(info):
    name = info2name(info[1:])
    alist = name.split('_')
    return alist

def info2name(info, addLabel=None):
    global badSp, unverifiedSp
    info = [each for each in info
            if each.strip() != '']
    sp = info[0]

    name = '_'.join(info)
    name = re.sub('\s', '_', name)

    if sp in badSp:
        name += '[unconfident_barcode!]'

    if sp in unverifiedSp:
        name = '[unverified]' + name

    if addLabel != None:
        if addLabel not in name:
            name += addLabel
    return name

def collape_identical_sequence(alist):
    if len(alist) <= 1:
        return alist

    blist = []
    taken_seq = set([])
    for name, seq in alist:
        if seq not in taken_seq:
            taken_seq.add(seq)
            blist.append((name, seq))

    #new: also take the one with fewer gap
    tmp = []
    bestGapN = 999
    lastLength = -1
    while(len(tmp) != lastLength):
        lastLength = len(tmp)
        tmp = []
        for ele in blist:
            name, seq = ele
            gapN = get_gapN(seq)
            if gapN <= bestGapN:
                tmp.append(ele)
                bestGapN = gapN

    blist = tmp
    return blist


def get_gapN(seq):
    try:
        gapN = gapN_dict[seq]

    except KeyError:
        gapN = seq.count('N') + seq.count('X') + seq.count('-')
        gapN_dict[seq] = gapN
    return gapN


def output_results(outdir, sorted_sps, dnlabel=''):
    ###the content in the main table
    ###This content include sp info and the barcode itself
    contents = []

    #the paired barcode with closest barcode found
    #TODO: make (assembled) disappear in itself but appear in pairs
    paired_fa = []

    #the barcode itself
    barcode_fa = []
    for sp in sorted_sps:
        info = read_spInfo(spInfo, sp)

        line = list(info)
        barcode = barcodes[sp]
        #TODO: also ordered by the label priority?
        barcode = collape_identical_sequence(barcode)
        line += [each[1] for each in barcode]
        contents.append(line)

        if len(barcode) > 1:
            dupSp.add(sp)
            if sp not in goodDup:
                print('warning: more than one barcode found for %s' % sp)

        for each in barcode:
            name, seq = each
            fasta = '>%s\n%s\n' % (info2name(read_spInfo(spInfo, sp), '(assembled)'), seq)
            paired_fa.append(fasta)
            barcode_fa.append(fasta)

        try:
            #could have more than one paired
            pairing = ''.join(['>%s\n%s\n' % (each[0], each[1]) for each in collape_identical_sequence(paired[sp])])
        except KeyError:
            #TODO: generate pairing for undo ones?
            pairing = '>NoPair_PlaceHolder\n%s\n' % ('-' * 658)
        paired_fa.append(pairing)

    print('unconfident samples: %s' % badSp)

    dn = '%s/%sverified_barcodes.fa' % (outdir, dnlabel)
    cmn.write_file(''.join(barcode_fa), dn)

    dn = '%s/%spaired_verified_barcodes.fa' % (outdir, dnlabel)
    cmn.write_file(''.join(paired_fa), dn)

    dn = '%s/%sall-barcodes.xlsx' % (outdir, dnlabel)
    workbook = xlsxwriter.Workbook(dn)
    worksheet = workbook.add_worksheet()
    #TODO: add format
    badSample = workbook.add_format({'bg_color': 'red'})
    dupSample = workbook.add_format({'bg_color': 'yellow'})
    unverSample = workbook.add_format({'bg_color': 'cyan'})
    notSureSample = workbook.add_format({'bg_color': 'blue'})
    dupNotSureSample = workbook.add_format({'bg_color': 'pink'})
    default = workbook.add_format()
    bold = workbook.add_format({'bold': True})

    worksheet.write(0, 0, '***unconfidnet barcode would be colored red', badSample)
    worksheet.write(1, 0, '***barcode with more than one form would be colored yellow', dupSample)
    worksheet.write(2, 0, '***barcode curated but still not sure would be colored blue', notSureSample)
    worksheet.write(3, 0, '***barcode requiring further manual curation would be colored cyan', unverSample)
    worksheet.write(4, 0, '***barcode of notSure and multiple copies would be colored pink', dupNotSureSample)

    headerI = 5 + 1
    for j, each in enumerate(header):
        worksheet.write(headerI, j, each, bold)
    worksheet.write(headerI, j+1, 'barcode', bold)


    rowStart = headerI + 1
    for i, line in enumerate(contents):
        row = rowStart + i
        sp = line[0]
        if sp in badSp:
            cellFormat = badSample
        elif (sp in notSureSet) and (sp in dupSp):
            cellFormat = dupNotSureSample
        elif sp in notSureSet:
            cellFormat = notSureSample
        elif sp in unverifiedSp:
            cellFormat = unverSample
        elif sp in dupSp:
            cellFormat = dupSample
        else:
            cellFormat = default
        for j, each in enumerate(line):
            worksheet.write(row, j, each, cellFormat)

    workbook.close()

def get_newest_tree_dir(alist):
    adict = {}
    for each in alist:
        print(each)
        name, date = each.split('_')[:2]
        dt = datetime.datetime.strptime(date, '%Y-%m-%d')
        value = (each, dt)
        try:
            adict[name].append(value)
        except KeyError:
            adict[name] = [value]

    rdict = {}
    for name in adict:
        alist = adict[name]
        if len(alist) == 1:
            rdict[name] = alist[0][0]
        else:
            bestOne = max(alist, key=lambda x: x[1])
            rdict[name] = bestOne[0]

    return rdict


def find_tree_file(wdir, name):
    tmp = cmn.cmd2lines('ls -t %s/*renamed' % wdir)
    if len(tmp) == 0:
        tmp += cmn.cmd2lines('ls -t %s/%s*tre' % (wdir, name))

    if len(tmp) == 0:
        print('Error! can not find tree in %s' % wdir)
        sys.exit()

    return tmp[0]


def order_ID_byTree(ftree, IDs):
    IDs = set(IDs)

    t = ete3.Tree(cmn.txt_read(ftree))

    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

    order_list = []
    for node in t:
        sp = node.name.split('_')[0].lstrip("'")
        if sp not in order_list and sp in IDs:
            order_list.append(sp)

    missing = set(IDs) - set(order_list)
    if len(missing) != 0:
        print('Error! missing IDs!')
        print('missed IDs: %s' % (', '.join(missing)))
        print('tree file: %s ' % ftree)
        sys.exit()
    return order_list



def split_and_order_sp_byTree(sampleIDs):
    treeDir = '/project/biophysics/Nick_lab/mtang/building_trees'
    wdirs = cmn.cmd2lines('ls %s' % treeDir)

    wdir_dict = get_newest_tree_dir(wdirs)
    rdict = {}
    for projectName in wdir_dict:
        wdir = '%s/%s' % (treeDir, wdir_dict[projectName])
        print('pick %s for %s' % (wdir, projectName))

        IDfile = '%s/NickList' % (wdir)
        IDlist = [each.replace('NVG-', '').replace('LEP-', 'LEP')
                for each in cmn.file2lines(IDfile)]
        overlapIDs = set(IDlist) & set(sampleIDs)
        if len(overlapIDs) > 0:
            ftree = find_tree_file(wdir, projectName)
            ordered_IDs = order_ID_byTree(ftree, overlapIDs)
            rdict[projectName] = ordered_IDs

    return rdict


def read_spInfo(spInfo, sp):
    global sp2name
    try:
        info = spInfo[sp]
    except KeyError:
        print('warnning: can not find info for %s, use info in the defline' % sp)
        info = sp2name[sp].split('_')
        if len(info) != 9:
            tmp = [info[0]]
            tmp += ['_'.join(info[1:])]
            tmp += [''] * 7
            info = tmp

    return info

def read_nickMade_barcodes(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue
        if line[0] == '#':
            continue

        if line[0] == '>':
            name = line[1:].strip().split('_')[0]

        seq = line.strip()
        adict[name] = seq
    return adict


if __name__=='__main__':

    fNick = '/project/biophysics/Nick_lab/wli/archive/barcodes/manually_curated_batches/Nick_curated_barcodes.txt'
    NickDict = read_nickMade_barcodes(fNick)
    NickSet = set(NickDict.keys())#used to record missing

    #wdir = '/project/biophysics/Nick_lab/wli/archive/barcodes'
    outdir = '/project/biophysics/Nick_lab/wli/archive/barcodes/auto_tables'

    fdup = '/project/biophysics/Nick_lab/wli/archive/barcodes/manually_curated_batches/verified_duplication.info'
    goodDup = set([line.strip().split()[0] for line in cmn.file2lines(fdup)])
    goodDup = goodDup | set(cmn.cmd2lines("grep assembled /project/biophysics/Nick_lab/wli/archive/barcodes/manually_curated_batches/unconfident_barcodes_paired.fa|cut -d '>' -f 2|cut -d '_' -f 1"))


    #spInfo is used in the final table
    #spNames is used as the defline for fasta and clw
    spInfo, header = gather_spInfo()

    #TODO: label unconfident barcode
    barcodeFiles = cmn.cmd2lines('ls -d /project/biophysics/Nick_lab/wli/archive/barcodes/manually_curated_batches/*| grep -v verified_duplication| grep -v old$|grep -v Readme| grep -v manual$|grep -v py$')
    print(barcodeFiles)
    #TODO: ad in the unconfident file, and also check for its existence in the normal barcodes to avoid mistakes
    #unconf_file = '/project/biophysics/Nick_lab/wli/archive/barcodes/manually_curated_batches/unconfident_barcodes_paired.fa'
    #barcodeFiles.remove(unconf_file)
    badIDs = set(cmn.file2lines('/project/biophysics/Nick_lab/wli/archive/barcodes/manually_curated_batches/unconfident_IDs.manual'))


    barcodes, paired = {}, {}
    dupSp = set([])
    #badSp = set([])
    badSp = badIDs - NickSet
    unverifiedSp = set([])
    notSureSet = set([])

    sp2name = {}#used to recover name for those with 'unknown'

    takenNickSet = set([])
    for fn in barcodeFiles:
        print('processing %s...' % fn)
        lines = [line for line in cmn.file2lines(fn)
                if line.strip() != '']

        #tell the type of file (clw or fasta)
        ordered_names, seqDict = readin_barcode_file(lines)

        #tell if it is paired
        #pairing could be done in two ways:
        #   1. one-to-one pairing with no label
        #   2. original reads are with 'assembled' and the others are paired
        pairType = tell_isPaired(ordered_names, fn)

        if pairType == False:
            print('guessing this file to be not paired')
        else:
            print('guessing this file to be paired (label: %s)' % pairType)

        sp = ''

        for i, name in enumerate(ordered_names):
            seq = seqDict[name]
            isPair = tell_pair_sequence(pairType, i, name)
            #print isPair, name
            fasta = (name, seq)
            if isPair: # is the paired sequence
                try:
                    paired[sp].append(fasta)
                except KeyError:
                    paired[sp] = [fasta]

            else: #original reads
                sp = name2sp(name)
                #print 'checkSp', sp, (sp in NickSet)
                if sp in NickSet:
                    takenNickSet.add(sp)
                    seq = NickDict[sp]
                    name = remove_labels(name)
                    fasta = (name, seq)

                sp2name[sp] = name
                if cmn.lastName(fn) == 'unconfident_barcodes_paired.fa' and (sp not in NickSet):
                    badSp.add(sp)

                if '[unverified]' in name and (sp not in NickSet):
                    unverifiedSp.add(sp)

                if '[notSure]' in name and (sp not in NickSet):
                    notSureSet.add(sp)
                try:
                    barcodes[sp].append(fasta)
                except KeyError:
                    barcodes[sp] = [fasta]


    #process output
    ###sort the output by the name to make genus together
    ###TODO: read the trees automatically to sort them
    sorted_sps = sorted(list(barcodes.keys()), key=lambda x: print4sort(read_spInfo(spInfo, x)))

    output_results(outdir, sorted_sps)

    #tree_groups = split_and_order_sp_byTree(sorted_sps)
    #subdir = '%s/barcodes_byTree' % outdir
    #cmn.mkdir(subdir)
    #for key in tree_groups:
    #    groupIDs = tree_groups[key]
    #    output_results(subdir, groupIDs, key+'_')

    #output finished IDs
    IDs = list(barcodes.keys())
    IDs = set(IDs) - badSp - unverifiedSp
    #print 'checkIDs'
    #print '\n'.join(IDs)
    lines = ['%s\t%s' % (ID, info2name(read_spInfo(spInfo, ID)))
            for ID in IDs]
    cmn.write_lines(lines, '%s/finished_sampleIDs.txt' % outdir)

    lines = ['%s\t%s' % (ID, info2name(read_spInfo(spInfo, ID)))
            for ID in badSp]
    cmn.write_lines(lines, '%s/unconfident_sampleIDs.txt' % outdir)

    lines = ['%s\t%s' % (ID, info2name(read_spInfo(spInfo, ID)))
            for ID in unverifiedSp]
    cmn.write_lines(lines, '%s/unverified_sampleIDs.txt' % outdir)

    cmd = 'cp %s %s' % (fNick, outdir)
    cmn.run(cmd)

    zipFile = 'allBarcodes.zip'
    cmd = 'rm %s/%s' % (outdir, zipFile)
    cmn.run(cmd)
    cmd = 'cd %s; zip -r %s *'  % (outdir, zipFile)
    cmn.run(cmd)

    #transfer
    cmd = 'rsync %s/* wenlin@morpho.swmed.edu:/data/www/wenlin/html/transfer/barcodes/' % outdir
    cmn.run(cmd)

    missing = NickSet - takenNickSet
    if len(missing) != 0:
        print('######################################################################')
        print('Important! the following samples Nick has looked but didn\'t get in:')
        print('\n'.join(missing))
        print('######################################################################')

    print('\nhttp://prodata.swmed.edu/wenlin/transfer/barcodes/allBarcodes.zip\n')

    #new, purge into the verification dataset
    print('purging barcodes into the verification dataset')
    fnew = '%s/verified_barcodes.fa' % outdir
    seqDict, length = read_fa(fnew)

    conn = bp.lock_database()

    madeOnes = bp.select_barcodes('isPCR=0', conn)
    for line in madeOnes:
        name = line[5]
        bp.delete_barcode_by_name(name, conn)
    conn.commit()


    for Oname in seqDict:
        name = Oname.replace('\'', '').replace('"', '').replace('(assembled)', '')
        name = name.split(']')[-1]
        print(name)
        seq = seqDict[Oname]
        sampleID = name.strip().split('_')[0]
        #sampleID, genus, sp = name.strip().split('_')[:3]
        #TODO: update names when sample info change?
        tmp, CSid = bp.find_IDs(name)
        bp.delete_barcode_by_name(name, conn, isCheck=False)
        bp.add_barcode(sampleID, 0, CSid, name, seq, conn)
    bp.close_database(conn)
