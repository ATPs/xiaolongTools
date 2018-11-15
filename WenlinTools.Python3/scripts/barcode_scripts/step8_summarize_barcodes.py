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
#import time
from collections import Counter
gapChars = set(['X', '-', 'N'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def check_conta_percent(fn):
    goodN = 0
    badN = 0
    count_good = True
    for line in cmn.file2lines(fn):
        if line[:4] == 'bait':
            continue
        if 'consensus' in line:
            continue
        if '############################' in line:
            count_good = False
            continue
        if count_good:
            goodN += 1
        else:
            badN += 1

    p = float(badN) / (badN + goodN)
    string = '%.2f(%s,%s)' % (p, badN, goodN)
    return string

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_inserted_gap(ID, seq, label):
    fn = 'sampleRun_%s/bait_insertion' % ID
    #if cmn.filexist(fn) or ('N' in seq.replace('-', 'N').strip('N')):
    if cmn.filexist(fn):
        #lines = cmn.file2lines(fn)
        #lines = sorted(lines, key=lambda x: int(x.split(',')[0][1:]))
        #Ngap = 0
        #for line in lines:
        #    items = line.strip().split()
        #    Ngap += len(items[-1])

        #check what is the right range of sequence
        print('runing blast to fix %s' % ID)
        checkSeq = seq.replace('-', 'N').strip('N')
        fquery = 'tmpInput.fa'
        fasta = '>input\n%s\n' % checkSeq
        cmn.write_file(fasta, fquery)
        dn = 'tmpBr_%s.txt' % label
        cmd = 'blastn -query %s -db /project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes_4verify.fa ' % fquery
        cmd += '-task blastn-short -dust no -outfmt \'6 qseqid sseqid qstart qend sstart send evalue pident qseq sseq\' -out %s' % dn
        cmn.run(cmd)
        isFixed = False
        for line in cmn.file2lines(dn):
            items = line.strip().split()
            #print items
            qstart, qend, sstart, send = list(map(int, items[2:6]))
            if sstart == 1 and send == 658 and qstart == 21:
                qseq, sseq = items[-2:]
                new = [char1 for char1, char2 in zip(qseq, sseq)
                        if char2 != '-']
                if len(new) == 658:
                    seq = seq[:qstart-1] + ''.join(new) + seq[qend:]
                    print('solution found for %s' % ID)
                    isFixed = True
                break
            if sstart == 2 and send == 655 and qstart == 22:
                qseq, sseq = items[-2:]
                new = [char1 for char1, char2 in zip(qseq, sseq)
                        if char2 != '-']
                if len(new) == 654:
                    seq = seq[:qstart-1] + ''.join(new) + seq[qend:]
                    print('solution found for %s' % ID)
                    isFixed = True
                break
        if not isFixed:
            cmn.append_file('%s\t%s\n' % (ID, label), 'cannot_fixed_indel.txt')
    return seq


def guessGenus(defline):
    subitems = defline.strip().split('_')
    genus = subitems[0]
    if any([char.isdigit() for char in genus]):
        genus = subitems[1]
    return genus

def isSameGenus(fn):
    lines = cmn.file2lines(fn)[1:]
    genus_list = []

    #There are two checks:
    #   1. if the denovo barcode has a same genus with ref, take it
    #   2. otherwise, check if denovo barcode has a same one as refbase
    #   3. if both not satisfied, then unknown
    for line in lines:
        items = line.strip().split()
        found_genus = guessGenus(items[3])
        ref_genus = guessGenus(items[4])
        N = int(items[2])
        if '_denovo' in items[0]:
            if found_genus == ref_genus:
                return ',takenD'
        if N <= 20:
            genus_list.append(found_genus)

    if len(genus_list) == len(lines):
        checkSet = set(genus_list)
        if len(checkSet) == 1:
            return ',takenD'
        elif len(checkSet) > 1:
            return ',diffGenus'
    return ',unknown'


def find_indel_from_reads(findel, fn):
    indel_list = get_indel_from_reads(fn)
    #indel_list = sum(indel_dict.values(), [])

    totalN = parse_indel_file(findel)
    del_p = []
    count_dict = Counter(indel_list)
    pStart_list = sorted(count_dict, key=lambda x: count_dict[x], reverse=True)

    takenN = 0
    for pStart, delN in pStart_list:
        for i in range(pStart, pStart + delN):
            del_p.append(i)

        takenN += delN
        if takenN == totalN:
            break
        elif takenN > totalN:
            print('Error! something went wrong in taking indel for %s' % fn)
            break

    #for delN in indel_inBaits:
    #    count = indel_inBaits[delN]
    #    pStart_list = find_indel_start(indel_dict, delN)
    #    for pStart in pStart_list[:count]:
    #        for i in xrange(pStart, pStart + delN):
    #            del_p.append(i)
    return set(del_p)

def find_indel_start(indel_dict, delN):
    alist = indel_dict[delN]
    count_dict = Counter(alist)
    sorted_list = sorted(count_dict, key=lambda x: count_dict[x], reverse=True)
    return sorted_list

def parse_indel_file(fn):
    #adict = {}
    total = 0
    for line in cmn.file2lines(fn):
        N = len(line.strip().split()[-1])
        total += N
        #try:
        #    adict[N] += 1
        #except KeyError:
        #    adict[N] = 1
    #return adict
    return total


def get_indel_from_reads(fn):
    #adict = {}
    alist = []
    seqLabels = ['b|', 'stack', 'thread']
    with open(fn) as fp:
        for line in fp:
            if any([label in line for label in seqLabels]):
                continue
            #find internal gaps
            seq = line.strip().split()[-1]
            isStart = False
            i = 0
            while(i < len(seq)):
                char = seq[i]
                if char.isupper():
                    isStart = True

                if char == '-' and isStart:
                    gapStart = i
                    gapLength = 1
                    alist.append((gapStart, 1))
                    while(i+1 < len(seq) and seq[i+1] == '-'):
                        i += 1
                        gapLength += 1
                        alist.append((gapStart - 1 + gapLength, 1))
                    #try:
                    #    adict[gapLength].append(gapStart)
                    #except KeyError:
                    #    adict[gapLength] = [gapStart]

                    #alist.append((gapStart, gapLength))
                i += 1
    return alist



def read_lineup_seq(fn, del_p):
    with open(fn) as fp:
        for line in fp:
            if 'stack' == line[:5]:
                stackSeq = line.strip().split()[-1]
            elif 'clean' == line[:5]:
                cleanSeq = line.strip().split()[-1]
            elif 'thread' == line[:6]:
                threadSeq = line.strip().split()[-1]

    stackSeq = parse_del_p(stackSeq, del_p)
    cleanSeq = parse_del_p(cleanSeq, del_p)
    threadSeq = parse_del_p(threadSeq, del_p)
    return threadSeq, stackSeq, cleanSeq

def parse_del_p(seq, del_p):
    new = [char for i, char in enumerate(seq)
            if i not in del_p]
    return ''.join(new)



if __name__=='__main__':
    #olines = cmn.cmd2lines('grep "thread\|stack" sampleRun_*/good_read_assembled.txt')
    #lines = cmn.cmd2lines('grep -H "threaded_\|stack_\|clean_" sampleRun_*/rescued_read_assembled_mis1*.txt')
    wdirs = cmn.cmd2lines('ls -d sampleRun_*')

    #cmn.mkdir('sampleRun_fake')
    #cmn.run('touch sampleRun_fake/good_read_assembled.txt')
    #cmn.run('touch sampleRun_fake/rescued_read_assembled_mis1.txt')
    #time.sleep(2)

    cmn.run('rm cannot_fixed_indel.txt 2> /dev/null')
    #frecords = cmn.cmd2lines('ls sampleRun_*/pickingLog.txt')

    stack_lines = {}
    thread_lines = {}
    clean_lines = {}
    leftN = 20
    barcodeLength = 658
    for wdir in wdirs:
        ID = wdir.split('sampleRun_')[-1]
        print('working on %s' % ID)
        try:
            fn = cmn.cmd2lines('ls sampleRun_%s/rescued_read_assembled_mis1*.txt' % ID)[0]
        except:
            print('can not find assembled files for %s' % ID)
            continue
        #print fn
        findel = 'sampleRun_%s/bait_insertion' % ID
        if cmn.filexist(findel):
            print('prasing indel for %s' % ID)
            indel_positions = find_indel_from_reads(findel, fn)
            print('indel_positions', indel_positions)
        else:
            indel_positions = []

        threadSeq, stackSeq, cleanSeq = read_lineup_seq(fn, indel_positions)

        thread_lines[ID] = threadSeq
        stack_lines[ID] = stackSeq
        clean_lines[ID] = cleanSeq

    #1. if thread and stack show inconsistent, show an X
    #2. if thread has gap, show as lower case
    #3. if both are gap, show an N
    #4. otherwise, show upper letter

    new = []
    refBaseDict = {}
    for Id in thread_lines:
        stackSeq = stack_lines[Id].replace('N', '-')[leftN: leftN+barcodeLength]
        threadSeq = thread_lines[Id].replace('N', '-')[leftN: leftN+barcodeLength]
        cleanSeq = clean_lines[Id].replace('N', '-')[leftN: leftN+barcodeLength]
        print(Id)
        #print stackSeq
        #print threadSeq

        seq = []
        for i, char1 in enumerate(cleanSeq):

            try:
                char2 = threadSeq[i]
            except:
                char2 = '-'

            try:
                char3 = stackSeq[i]
            except:
                char3 = '-'

            if char1 == char2:
                seq.append(char2)
            else:#different characters
                if char1 == '-' and char2 == '-':
                    seq.append(char3.lower())
                elif char1 == '-':
                    seq.append(char2.lower())
                elif char2 == '-':
                    seq.append(char1.lower())
                else:
                    #different chars and not a gap
                    seq.append('X')
        fasta = '>%s\n%s\n' % (Id, ''.join(seq))
        refBaseDict[Id] = ''.join(seq)
        new.append(fasta)

    cmn.write_file(''.join(new), 'sum_barcodes.fa')

    #cmn.run('rm -r sampleRun_fake')

    #check denovo pipeline one
    fns = cmn.cmd2lines('ls sampleRun_*/denovo_barcode.fa')
    denovoDict = {}
    new = []
    for fn in fns:
        Id = cmn.find_between(fn, 'sampleRun_', '/')
        lines = cmn.file2lines(fn)
        seq = ''.join(lines[1:])
        if seq > 658:
            tmp = seq.replace('N', '')
            if len(tmp) == 658:
                seq = tmp
        denovoDict[Id] = seq

        fasta = '>%s\n%s\n' % (Id, seq)
        new.append(fasta)

    cmn.write_file(''.join(new), 'sum_denovo.fa')

    new = []
    for Id in clean_lines:
        seq = clean_lines[Id]
        fasta = '>%s\n%s\n' % (Id, seq[20:678])
        new.append(fasta)

    cmn.write_file(''.join(new), 'sum_clean.fa')


    #checking concamination
    #fns = cmn.cmd2lines('ls sampleRun_*/good_read_assembled.txt')
    fns = cmn.cmd2lines('ls sampleRun_*/rescued_read_assembled_mis1*.txt')
    conta_dict = {}
    for fn in fns:
        Id = cmn.find_between(fn, 'sampleRun_', '/')
        percent = check_conta_percent(fn)
        conta_dict[Id] = percent

    #compare barcodes
    IDs = [each.split('_')[-1] for each in wdirs]

    #15099E08        hasRef  Gap0    completeDenovo  diff_2  0.18(229,1015)
    header = 'sample refInfo refGap denovoInfo diffLabel conta_fration(conta:good)'
    new = ['\t'.join(header.strip().split())]
    for ID in IDs:
        line = [ID]
        frecord = 'sampleRun_%s/pickingLog.txt' % ID
        recordLabel = ''
        if cmn.filexist(frecord):
            record_info = cmn.txt_read(frecord)
            if 'noGenus' in record_info:
                recordLabel = 'noGenus'
            elif 'noSpecies' in record_info:
                recordLabel = 'noSpecies'

        try:
            refSeq = refBaseDict[ID].upper()
            line.append('hasRef')
        except KeyError:
            refSeq = None
            line.append('noRef')

        if refSeq == None:
            line.append('NA')
        else:
            N = sum([refSeq.count(char) for char in gapChars])
            line.append('Gap%s' % N)

        try:
            denovoSeq = denovoDict[ID].upper()
            if 'N' in denovoSeq or ('-' in denovoSeq) or len(denovoSeq) < 658:
                line.append('inCompleteDenovo%s' % len(denovoSeq.replace('N','')))
            else:
                tmplabel = 'completeDenovo'
                try:
                    fcheck = 'sampleRun_%s/%s_summaryCheck.table' % (ID, ID)
                    sublabel = isSameGenus(fcheck)
                    if sublabel == ',unknown':
                        N = sum([refSeq[i] != denovoSeq[i]
                            for i in range(len(refSeq))
                            if denovoSeq[i] not in gapChars])
                        if N <= 10:
                            sublabel = ',takenD'
                    tmplabel += sublabel
                except IOError:
                    pass
                line.append(tmplabel)
        except KeyError:
            denovoSeq = None
            line.append('noDenovo')

        #construct the label
        difflabel = ''
        if refSeq == None:
            difflabel = 'Error'
        elif denovoSeq == None:
            difflabel = 'hardCase'
        else:
            if refSeq == denovoSeq:
                difflabel = 'same,takenD'
            else:
                if len(refSeq) != len(denovoSeq):
                    difflabel = 'diffLength'
                    #difference = compare_sequence(refSeq, denovoSeq)
                    #difflabel += 'D%s' % difference
                else:
                    N = sum([refSeq[i] != denovoSeq[i]
                        for i in range(len(refSeq))
                        if denovoSeq[i] not in gapChars])

                    difflabel = 'diff_%s' % N
                    #check genus, if they are in the same genus, suggest taken

        try:
            clean_seq = clean_lines[ID]
            thread_seq = thread_lines[ID]
            length = min(len(clean_seq), len(thread_seq))
            N = sum([clean_seq[i+20] != thread_seq[i+20] for i in range(length-20)
                if clean_seq[i+20] not in gapChars and (thread_seq[i+20] not in gapChars)])

            Ngap = sum([char in gapChars for char in clean_seq[20:678]])
            if Ngap == 0 and N == 0:
                difflabel += ',goodCC'
            else:
                difflabel += ',CC_%s' % N
                if Ngap != 0:
                    difflabel += '[g%s]' % Ngap
        except KeyError:
            difflabel += ',noCC'

        line.append(difflabel)
        try:
            line.append(conta_dict[ID])
        except:
            line.append('NA')
        line.append(recordLabel)
        new.append('\t'.join(line))
    new.append('')
    cmn.write_lines(new, 'compare.check')
