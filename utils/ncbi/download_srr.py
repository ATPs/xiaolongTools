# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:34:58 2019
download srr
@author: ATPs
"""
import os

def srr2ftploc(srr):
    '''
    given a srr, return the location in ncbi ftp
    '''
    return 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'+srr[:3]+'/'+srr[:6]+'/'+srr+'/'+srr+'.sra'

def getSRRip(srr, workfolder = '.'):
    '''
    input a srr, download srr.sra in workfolder
    default workfolder is current one
    '''
    cmd = 'cd '+workfolder +' && wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'+srr[:3]+'/'+srr[:6]+'/'+srr+'/'+srr+'.sra'
    os.system(cmd)

description = '''input a srr, download srr.sra in workfolder
    default workfolder is current one
    download_srr.py -i SRR1576990 -o /home/x/RNA
    
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file fasta file', required=True)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored',default=None)
    f = parser.parse_args()
    getSRRip(srr=f.input, workfolder=f.output)

