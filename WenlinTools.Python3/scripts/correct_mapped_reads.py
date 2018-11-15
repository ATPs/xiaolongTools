#!/home2/wli/anaconda/bin/python

#function:
#input:
#output:
#algorithm: use the mapping relationship to correct SNP and make it to the reference
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tdict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
        }

def reverse_strand(seq):
    new = []
    for char in seq[::-1]:
        try:
            new.append(tdict[char])
        except KeyError:
            new.append('N')
    return ''.join(new)


def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0].split()[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def correct_mapSeq(seq, refSeq, pairing):
    #cmpSeq = seq[qstart:qend]
    new = []
    for pair in pairing:
        i, refI = pair
        if i == None:
            continue

        if refI == None:
            new.append(seq[i])
        else:
            new.append(refSeq[refI])

    new = ''.join(new)

    return new


def recover_fq(mRead, newSeq):
    if mRead.is_read1:
        label = 'R1'
    elif mRead.is_read2:
        label = 'R2'
    else:
        label = 'singleton'

    qual = mRead.qual
    if mRead.is_reverse:
        newSeq = reverse_strand(newSeq)
        qual = qual[::-1]

    ID = '@%s' % mRead.query_name
    if label == 'R1':
        ID += '/1'
    elif label == 'R2':
        ID += '/2'
    new = [ID]
    new.append(newSeq)
    new.append('+')
    new.append(qual)
    new.append('')

    newfq = '\n'.join(new)
    return newfq, label



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fref = sys.argv[1:3]
    except:
        print("Usage: *.py *.sam ref.fa", file=sys.stderr)
        sys.exit()

    seqDict = read_fa(fref)

    samfile = pysam.AlignmentFile(fn)

    rdict = {}
    dp_unmapped = open(cmn.lastName(fn).replace('.sam', '') + '_unmapped.fastq', 'w')

    for mRead in samfile:
        if mRead.is_secondary:
            continue

        if mRead.is_unmapped:
            fq, dlabel = recover_fq(mRead, mRead.query_sequence)
            dp_unmapped.write(fq)
            continue

        scaf = mRead.reference_name
        refSeq = seqDict[scaf]
        #refHsp = refSeq[mRead.reference_start: mRead.reference_end]
        newSeq = correct_mapSeq(mRead.seq, refSeq, mRead.aligned_pairs)
        #print mRead.query_name
        #print 'refHsp', refHsp
        #print 'newSeq', newSeq
        fq, dlabel = recover_fq(mRead, newSeq)
        try:
            rdict[dlabel].append(fq)
        except:
            rdict[dlabel] = [fq]

    samfile.close()
    dp_unmapped.close()

    sp = cmn.lastName(fn).split('_')[0]
    for dlabel in rdict:
        dn = '%s_%s_corrected.fastq' % (sp, dlabel)
        #cmd = 'grep ^@ %s > %s' % (fn, dn)
        #cmn.run(cmd)
        dp = open(dn, 'w')
        for each in rdict[dlabel]:
            dp.write(each)
        dp.close()


