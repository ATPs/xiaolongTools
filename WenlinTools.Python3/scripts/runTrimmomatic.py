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


reverse_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        }

def reverse_strand(read):
    new = []
    for char in read[::-1]:
        try:
            a = reverse_dict[char]
        except KeyError:
            a = 'N'
        new.append(a)

    return ''.join(new)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def output_fastq(dnR1, dnR2, fq1, fq2):
    for dn, fq in [(dnR1, fq1), (dnR2, fq2)]:
    #for dn, fq in [(dnR1, fq1)]:
        for line in fq:
            dn.write(line)


def output_fastq_merge(dn, fq):
    #for dn, fq in [(dnR1, fq1), (dnR2, fq2)]:
    #for dn, fq in [(dnR1, fq1)]:
    if True:
        for line in fq:
            dn.write(line)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def check_trim_log(fn):
    aset, bset = set([]), set([])

    unpaired = set([])
    fp = open(fn)
    for line in fp:
        #J00138:107:HJMTYBBXX:2:1101:17624:1015 1:N:0:NACCTGC 0 0 0 0
        #J00138:107:HJMTYBBXX:2:1101:17624:1015/1
        items = line.strip().split()
        ID = items[0]
        fqLabel = None
        if '/' in ID:
            a, b = ID.split('/')
            ID = a
            fqLabel = b
        if ID in unpaired:
            continue
        slength = int(items[-4])
        if slength == 0:
            unpaired.add(slength)
        else:
            trimN = int(items[-1])
            if trimN != 0:
                if items[1][0] == '1' or fqLabel == '1':
                    aset.add(ID)
                elif items[1][0] == '2' or fqLabel == '2':
                    bset.add(ID)
                else:
                    print('Error: something wrong with %s' % line)
    fp.close()
    return aset, bset


def find_match_part(seq1, seq2):
    #need to reverse one ot them
    minlength = min([len(seq1), len(seq2)])
    #if minlength < 30:
    #    seed_length = minlength - 3
    #else:
    #    seed_length = minlength - 10
    seed_length = 20
    seed_length = min(minlength, 20)
    rseq1 = reverse_strand(seq1)
    i1 = -1
    i2 = -1
    i = 0

    #find seed match
    #print 'pre', i, seed_length, len(rseq1)
    while((i + seed_length) <= len(rseq1)):
        checkSeq = rseq1[i:i+seed_length]
        tmpI = seq2.find(checkSeq)
        #print 'interM', i, tmpI
        if tmpI != -1:
            i1 = i
            i2 = tmpI
            break
        i += 1

    if tmpI == -1:
        print('trimInfo:', ID, 'noMatch')
        return None, None, None, None

    j1 = i1
    j2 = i2
    #extend to the end
    while(i1 > 0 and i2 > 0):
        i1 -= 1
        i2 -= 1

    while(j1 < len(rseq1) and j2 < len(seq2)):
        j1 += 1
        j2 += 1

    #reverse the index for seq1
    i1 = len(rseq1) - i1
    j1 = len(rseq1) - j1
    i1, j1 = j1, i1
    print('trimInfo:', ID, i, tmpI, len(seq1), i1, j1, len(seq2), i2, j2)
    return i1, j1, i2, j2


def trim_fq(i, j, fq):
    defline, seq, sth, qc = fq
    new = [defline, seq[i:j].strip() + '\n', sth, qc[i:j].strip() + '\n']
    return new

if __name__=='__main__':
    #options=parse_options()
    try:
        fR1, fR2, ad1, outlabel = sys.argv[1:]
    except:
        print("Usage: *.py R1 R2 ad1 outlabel", file=sys.stderr)
        sys.exit()


    #prepare adaptor
    dn_ad = './%s_adapter.fa' % outlabel
    fasta = '>ad/1\n%s\n' % ad1
    fasta += '>ad/2\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n'

    cmn.write_file(fasta, dn_ad)

    cmd = 'java -Xmx10g -jar /home2/wli/local/Trimmomatic-0.36/trimmomatic-0.36.jar PE %s %s ' % (fR1, fR2)
    cmd += '-baseout %s ILLUMINACLIP:%s:4:10:10:8:true ' % (outlabel, dn_ad)
    cmd += ' -trimlog %strim.log ' % outlabel
    print('running command: %s' % cmd)
    cmn.run(cmd)


    #here, apply second filter to ensure the paired ones are reverse complement
    #if not, trim out the extended part
    print('begin post processing...')
    isTrimR1, isTrimR2 = check_trim_log(outlabel + 'trim.log')
    print('finish reading trimming info...')

    fnR1 = '%s_1P' % outlabel
    fnR2 = '%s_2P' % outlabel
    dnR1 = '%shc_R1.fastq' % outlabel
    dnR2 = '%shc_R2.fastq' % outlabel
    dnMerge = '%shc_merge.fastq' % outlabel

    NnoMatch = 0
    Ntotal = 0
    Nprocess = 0
    with open(fnR1) as fp1, open(fnR2) as fp2, open(dnR1, 'w') as dn1, open(dnR2, 'w') as dn2, open(dnMerge, 'w') as dm:
        for i, line1 in enumerate(fp1):
            line2 = fp2.readline()
            if i % 4 == 0:
                ID = line1[1:].split()[0].split('/')[0]
                if (ID in isTrimR1) or (ID in isTrimR2):
                    isNeedCheck = True
                else:
                    isNeedCheck = False
                fq1 = [line1]
                fq2 = [line2]

            else:
                fq1.append(line1)
                fq2.append(line2)

            if i % 4 == 3:
                Ntotal += 1
                #single output

                #if isNeedCheck:
                if True:
                    seq1 = fq1[1].strip()
                    seq2 = fq2[1].strip()
                    #if len(seq1) == len(seq2):
                    if False:
                        #good one, just output
                        output_fastq(dn1, dn2, fq1, fq2)
                    else:
                        Nprocess += 1
                        i1, j1, i2, j2 = find_match_part(seq1, seq2)
                        if i1 == None:
                            print('skip nonMatch read pair: %s' % ID)
                            NnoMatch += 1
                            continue

                        if j1 - i1 < len(seq1) * 0.8 and j2 - j2 < len(seq2) * 0.8:
                            output_fastq(dn1, dn2, fq1, fq2)
                            continue

                        fq1 = trim_fq(i1, j1, fq1)
                        fq1[0] = '@%s\n' % ID
                        #fq2 = trim_fq(i2, j2, fq2)
                        output_fastq_merge(dm, fq1)

                else:
                    output_fastq(dn1, dn2, fq1, fq2)
                    if abs(len(fq1[1]) - len(fq2[1])) > 20:
                        print('Warnning: this read length differ too much: %s' % ID)

    info = 'total: %s\nprocess: %s\nnoMatch: %s\n' % (Ntotal, Nprocess, NnoMatch)
    cmn.write_file(info, outlabel + '.statInfo')
