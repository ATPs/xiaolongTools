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
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, ref_scaf = sys.argv[1:]
    except:
        print("Usage: *.py sam scaffold", file=sys.stderr)
        sys.exit()


    samfile = pysam.AlignmentFile(fn)

    #limit = [0, 80]

    aln_list = set([])
    nameDict = {}
    for record in samfile:
        if not record.is_unmapped and (not record.is_secondary):
            scaf = record.reference_name
            if scaf != ref_scaf:
                continue
            seq = record.seq
            aligned = record.aligned_pairs
            subjct_range = [i[1] for i in aligned if i[1]!=None]
            left = min(subjct_range)
            right = max(subjct_range)
            #change unmapped letter into lower case

            seq = list(seq)
            for i, j in aligned:
                if j == None:
                    seq[i] = seq[i].lower()
            seq = ''.join(seq)


            #if left > limit[1] or right < limit[0]:
            #    continue

            #print record.query_name,left, right, aligned
            if aligned[0][1] == None:
                #need to fill in beginning
                shifts = []
                for i, each in enumerate(aligned):
                    if each[1] == None:
                        shifts.insert(0, -i - 1)
                    else:
                        break

                print('before', left, aligned)
                for i, shift in enumerate(shifts):
                    aligned[i] = (aligned[i][0], left + shift)
                left -= len(shifts)
                print('after', left, aligned)

                if left < 0:#shift too much
                    for i in range(0 - left):
                        aligned[i] = (aligned[i][0], None)

                    left = 0
                    print('fix', left, aligned)


            aln = []
            for _ in range(left):
                aln.append('-')

            isPass = False
            for i, j in aligned:
                if j == None:
                    if isPass:
                        aln.append(seq[i])
                    continue
                if i == None:
                    aln.append('-')
                    continue
                isPass = True
                aln.append(seq[i])

            aln = ''.join(aln)
            aln_list.add(aln)
            #aln_list.add(aln)
            name = ':'.join(record.query_name.split(':')[-3:])
            if record.is_read1:
                name += '_1'
            elif record.is_read2:
                name += '_2'

            try:
                nameDict[aln].append(name)
            except KeyError:
                nameDict[aln] = [name]

    dn = cmn.lastName(fn).replace('.sam', '') + '_aln.txt'
    info = []
    maxLength = max([len(each) for each in aln_list])
    maxNameLength = max([len('|'.join(each)) for each in list(nameDict.values())])
    nameformat = '{:<%s}' % maxNameLength
    for i, aln in enumerate(aln_list):
        #name = 'readgroup%s' % i
        name = nameformat.format('|'.join(nameDict[aln]))

        toAdd = maxLength - len(aln)
        if toAdd > 0:
            aln += '-' * toAdd
        info.append('%s    %s\n' % (name, ''.join(aln)))
    cmn.write_file(''.join(info), dn)
