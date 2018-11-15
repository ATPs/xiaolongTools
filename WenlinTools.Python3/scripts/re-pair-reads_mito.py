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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def separate_reads(fn):
    adict = {}
    bdict = {}

    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 0:
                ID = line[1:].strip()
                record = []

            record.append(line)

            if i % 4 == 3:
                key = ID.split('/')[0]
                if ID[-1] == '1':
                    adict[key] = ''.join(record)
                elif ID[-1] == '2':
                    bdict[key] = ''.join(record)

    return adict, bdict




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, outlabel = sys.argv[1:3]
    except:
        print("Usage: *.py in.fastq outlabel", file=sys.stderr)
        sys.exit()


    reads1, reads2 = separate_reads(fn)

    read1_keys = set(reads1.keys())
    read2_keys = set(reads2.keys())

    paired_keys = read1_keys & read2_keys

    r1, r2, single = [], [], []

    for key in paired_keys:
        fq1 = reads1[key]
        r1.append(fq1)
        fq2 = reads2[key]
        r2.append(fq2)

    unpaired_R1 = read1_keys - paired_keys
    for key in unpaired_R1:
        fq = reads1[key]
        single.append(fq)

    unpaired_R2 = read2_keys - paired_keys
    for key in unpaired_R2:
        fq = reads2[key]
        single.append(fq)


    dn = outlabel + '_paired1.fastq'
    cmn.write_file(''.join(r1), dn)

    dn = outlabel + '_paired2.fastq'
    cmn.write_file(''.join(r2), dn)

    dn = outlabel + '_unpaired.fastq'
    cmn.write_file(''.join(single), dn)


