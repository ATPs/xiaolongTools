#!/usr6/local/bin/python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import cmn
from collections import Counter
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fastq(fn):
    adict = {}
    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 0:
                ID = line.strip()
                record = [line]
                continue

            #other lines would be appended
            record.append(line)

            if i % 4 == 3:
                adict[ID] = record
    return adict


def combine_fastq(oldDict, newDict):
    global ischeck_badID
    IDs = set(oldDict.keys()) | set(newDict.keys())
    #1. check the first part of ID to see if they are the same
    if ischeck_badID:
        checkIDs = [ID.split()[0] for ID in IDs]
        count_dict = Counter(checkIDs)
        bad_IDs = [ID for ID in checkIDs if count_dict[ID] != 1]
        if len(bad_IDs) != 0:
            print('Error! duplicated machine ID found with different affix!')
            print('here are the bad machine IDs:')
            print('\n'.join(bad_IDs))
            print('\nIf you believe these are good, you can append \'False\' in your command to suppress this check')
            sys.exit()
    rdict = {}
    for ID in IDs:
        records = []
        try:
            oldR = oldDict[ID]
            records.append(oldR)
        except KeyError:
            pass

        try:
            newR = newDict[ID]
            records.append(newR)
        except KeyError:
            pass

        if len(records) == 1: # no duplication
            rdict[ID] = ''.join(records[0])
        else:#has duplicaiton
            r1, r2 = records
            if len(r1[1]) > len(r2[1]):
                r = r1
            else:
                r = r2
            rdict[ID] = ''.join(r)
    return rdict


def combine_data(fn, label):
    global ischeck_badID, sp, old_dir
    sp = '_'.join(cmn.lastName(fn).split('_')[:-1])

    accepted_labels = ['R1', 'R2', 'singleton']
    if label not in accepted_labels:
        print('Error! your indicated label are not accepted')
        print('accepted values are %s' % (','.join(accepted_labels)))
        sys.exit()

    newDict = read_fastq(fn)

    oldFn = '%s/%s_%s.fastq' % (old_dir, sp, label)
    if cmn.filexist(oldFn):
        print('combine new data with old data for %s' % fn)
        oldDict = read_fastq(oldFn)
        #newDict = read_fastq(fn)
        finalDict = combine_fastq(oldDict, newDict)
    else:
        finalDict = newDict

    return finalDict


def repair_reads(adict, bdict):
    print('pairing reads...')
    key1s = {each.split('/')[0].split()[0]: each for each in adict}
    key2s = {each.split('/')[0].split()[0]: each for each in bdict}
    keys = set(key1s.keys()) & set(key2s.keys())

    if len(keys) != len(key1s) or (len(keys) != len(key2s)):
	    print('Error! R1 and R2 has unpaired reads! Please check!')
	    sys.exit()
    keys = list(keys)
    alist, blist = [], []
    for key in keys:
        record1 = adict[key1s[key]]
        alist.append(record1)
        record2 = bdict[key2s[key]]
        blist.append(record2)
    return alist, blist

def check_fastqlines(fn):
    cmd = 'wc -l %s' % fn
    N = int(cmn.cmd2info(cmd).strip().split()[0])
    return N

if __name__=='__main__':
    #options=parse_options()
    try:
        fn1, fn2 = sys.argv[1:3]
    except:
        print("Usage: *.py R1.fastq R2.fastq", file=sys.stderr)
        sys.exit()


    old_dir = '/project/biophysics/Nick_lab/mtang/archive/fastq_libs'

    sp1 = '_'.join(cmn.lastName(fn1).split('_')[:-1])
    sp2 = '_'.join(cmn.lastName(fn2).split('_')[:-1])
    if sp1 != sp2:
        print('Error! the sample ID is not the same for R1 and R2!')
        sys.exit()
    else:
        sp = sp1
    
    #if len(cmn.lastName(fn1).split('_')) > 2:
    #    print 'Error! the code is not designed to handle such libs'
    #    print 'Please update it manually'
    #    sys.exit()


    ischeck_badID = True
    try:
        clabel = sys.argv[3]
        if clabel == 'False':
            print('suppress machine ID check')
            ischeck_badID = False
    except:
        pass

    oldFn1 = '%s/%s_R1.fastq' % (old_dir, sp)
    oldFn2 = '%s/%s_R2.fastq' % (old_dir, sp)

    if (not os.path.exists(oldFn1)) and (not os.path.exists(oldFn2)):
        print('fastq is not present in destination directory, copy them in directly')
        cmd = 'cp %s %s' % (fn1, oldFn1)
        cmn.run(cmd)
        print('your fastq is updated into %s' % oldFn1)
        cmd = 'cp %s %s' % (fn2, oldFn2)
        cmn.run(cmd)
        print('your fastq is updated into %s' % oldFn2)
        sys.exit()
    else:
        #check if they have the same length of lines
        #if so, thinking they are the same
        N1 = check_fastqlines(fn1)
        N1o = check_fastqlines(oldFn1)
        N2 = check_fastqlines(fn2)
        N2o = check_fastqlines(oldFn2)
        if N1 == N1o and N2 == N2o:
            print('the file %s is the same as %s, skip' % (fn1, oldFn1))
            print('the file %s is the same as %s, skip' % (fn2, oldFn2))
            sys.exit()


    fn1_dict = combine_data(fn1, 'R1')
    fn2_dict = combine_data(fn2, 'R2')

    fn1_list, fn2_list = repair_reads(fn1_dict, fn2_dict)

    oldFn = '%s/%s_R1.fastq' % (old_dir, sp)
    info = ''.join(fn1_list)
    cmn.write_file(info, oldFn)
    print('your fastq is updated into %s' % oldFn)

    oldFn = '%s/%s_R2.fastq' % (old_dir, sp)
    info = ''.join(fn2_list)
    cmn.write_file(info, oldFn)
    print('your fastq is updated into %s' % oldFn)



