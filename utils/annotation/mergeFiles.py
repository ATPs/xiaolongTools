# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 11:02:07 2019

@author: ATPs
"""

def mergefiles(files_in, file_out):
    '''
    files_in is files to merge, file_out is where to store the output file
    files_in is a string that will be used for glob.glob match
    '''
    import glob
    files = glob.glob(files_in)
    print('totally ',len(files),'files to merge')
    if len(files)>1:
        print('the first and last files are', files[0],files[-1])
    fout = open(file_out,'w')
    for f in files:
        fout.write(open(f).read())
    fout.close()

description = '''
    files_in is files to merge, file_out is where to store the output file
    files_in is a string that will be used for glob.glob match
'''
if __name__ == '__main__':
    print(description)
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--files_in', help = 'a string that will be used for glob.glob match', required=True)
    parser.add_argument('-o','--file_out', help = 'file_out is where to store the output file, default merged.result', default='merged.result')
    f = parser.parse_args()
    mergefiles(files_in=f.files_in, file_out=f.file_out)