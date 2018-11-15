#!/usr/bin/env python3

import os

def splitFile2Nparts(filename, N = 10, outfolder = None):
    '''
    Given a filename of file with many lines, split the file to N files
    each file with almost equal size
    if outfolder is given, then put outfiles in that folder.
    the output file name will be filename+'.split0'
    '''
    filefolder = os.path.dirname(filename)
    filebasename = os.path.basename(filename)
    
    #change to working directory
    if outfolder is not None:
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
            os.chdir(outfolder)
    else:
        os.chdir(filefolder)
    
    filesOut = [open(filebasename+'.split'+str(i), 'w') for i in range(N)]
    for n, line in enumerate(open(filename,'r')):
        filen = n % N
        filesOut[filen].write(line)
    for f in filesOut:
        f.close()
    print('done')

description = '''
Given a filename of file with many lines, split the file to N files
each file with almost equal size
if outfolder is given, then put outfiles in that folder.
the output file name will be filename+'.split0'
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input the file that need to be splitted', required=True)
    parser.add_argument('-N','--Number',help = 'number of parts to split', default = 10,type=int)
    parser.add_argument('-o','--output',help = 'folder of the where the output files stored', default = None)
    f = parser.parse_args()
    splitFile2Nparts(filename = f.input, N = f.Number, outfolder = f.output)
