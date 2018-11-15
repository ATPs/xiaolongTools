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
        print("Usage: *.py fasta_file", file=sys.stderr)
        sys.exit()


    #f_finished = '/project/biophysics/Nick_lab/wli/sequencing/verify_barcodes/finished_barcodes'
    f_finished = '/project/biophysics/Nick_lab/wli/sequencing/verify_barcodes/auto_tables/finished_sampleIDs.txt'

    finished_IDs = set([each.split()[0].replace('[new]', '') for each in cmn.file2lines(f_finished)])

    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                line = line.strip()
                spID = line[1:].split('_')[0]
                if spID in finished_IDs:
                    print('finished:', line)




