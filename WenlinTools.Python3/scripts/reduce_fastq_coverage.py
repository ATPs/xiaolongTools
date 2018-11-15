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
        fR1, fR2 = sys.argv[1:3]
    except:
        print("Usage: *.py R1 R2", file=sys.stderr)
        sys.exit()


    sample = cmn.lastName(fR1).split('_')[0]

    spacing_list = [2, 3, 5, 10]

    count = 0
    for spacing in spacing_list:
        fpR1 = open(fR1)
        fpR2 = open(fR2)
        dnlabel = '%s.spacing%s' % (sample, spacing)
        print('making %s' % dnlabel)

        dnR1 = open('%s_R1.fastq' % dnlabel, 'w')
        dnR2 = open('%s_R2.fastq' % dnlabel, 'w')

        for i, line1 in enumerate(fpR1):
            line2 = fpR2.readline()
            if i % 4 == 0:
                count += 1

            if count % spacing == 1:
                dnR1.write(line1)
                dnR2.write(line2)

        dnR1.close()
        dnR2.close()
        fpR1.close()
        fpR2.close()



