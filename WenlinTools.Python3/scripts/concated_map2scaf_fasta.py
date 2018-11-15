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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def output_result(scaf, result):
    dn = scaf + '_cm2f.fa'
    with open(dn, 'w') as dp:
        for name in result:
            fasta = '>%s\n%s\n' % (name, ''.join(result[name]))
            dp.write(fasta)





if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py concat.map (this map should have header)", file=sys.stderr)
        sys.exit()


    result = {}
    lastScaf = None
    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i == 0:
                header = line.strip().split()[2:]
            else:
                items = line.strip().split()
                scaf = items[0]

                if lastScaf == None:
                    lastScaf = scaf
                else:
                    if lastScaf != scaf:
                        #output result
                        output_result(lastScaf, result)
                        result = {}
                        lastScaf = scaf

                for j, item in enumerate(items[2:]):
                    name = header[j]
                    try:
                        result[name].append(item)
                    except KeyError:
                        result[name] = [item]

    output_result(lastScaf, result)

