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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py fastq", file=sys.stderr)
        print("Assuming the +33 scale", file=sys.stderr)
        sys.exit()

    dn = fn + '.corrected'
    with open(dn, 'w') as dp:
        with open(fn) as fp:
            for i, line in enumerate(fp):
                if i % 4 == 3:
                    qletters = line.strip()
                    Min = ord(min(qletters))
                    Max = ord(max(qletters))
                    if Max > 74:
                        #assuming the +64 scale
                        if Min < 64:
                            print('Error! weird line found: max %s, min %s' % (Max, Min))
                            print('line: %s' % qletters)
                        else:
                            new = [chr(ord(each) - 31) for each in qletters]
                            print('corrected one line, max %s, min %s' % (Max, Min))
                            dp.write(''.join(new))
                            dp.write('\n')
                    else:
                        #normal +33 scale
                        dp.write(line)
                else:
                    dp.write(line)



