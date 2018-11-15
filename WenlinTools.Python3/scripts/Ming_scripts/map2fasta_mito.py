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
transfer_dict = {
        'CT': 'R',
        'AG': 'Y',
        'AT': 'W',
        'CG': 'S',
        'GT': 'M',
        'AC': 'K'
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.map", file=sys.stderr)
        sys.exit()


    label = cmn.lastName(fn).split('_')[0]
    dn = cmn.lastName(fn).replace('.map', '') + '_m2s.fa'

    print('parsing for %s' % label)
    print('would save result to %s' % dn)

    a = []
    b = []

    with open(fn) as fp:
        for line in fp:
            try:
                i, j = line.strip().split()
            except:
                i = line.strip()
                j = None
            a.append(i)
            b.append(j)

    if j == None:#only one line in map
        fasta = '>%s_mito\n%s\n' % (label, ''.join(a))
    else:
        #check if the two sequence are the same
        diff_dict = {}
        for i, aa in enumerate(a):
            bb = b[i]
            if aa != bb:
                diff_dict[i] = [aa, bb]
        if len(diff_dict) == 0:
            fasta = '>%s_mito\n%s\n' % (label, ''.join(a))
        else:
            print('polymorphism detected for the following positions:')
            for i in diff_dict:
                print(i+1, diff_dict[i])
            print('will use ambiguous bp for those positions')
            
            ambP = set(diff_dict.keys())
            seq = []
            for i, aa in enumerate(a):
                if i not in ambP:
                    seq.append(aa)
                else:
                    bps = diff_dict[i]
                    bps.sort()
                    code = transfer_dict[''.join(bps)]
                    seq.append(code)
            fasta = '>%s_mito\n%s\n' % (label, ''.join(seq))        
                    
        #fasta = '>%s_cp1\n%s\n>%s_cp2\n%s\n' % (label, ''.join(a), label, ''.join(b))
    cmn.write_file(fasta, dn)



