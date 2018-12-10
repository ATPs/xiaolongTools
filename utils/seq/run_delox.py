# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 20:56:01 2018

@author: ATPs
"""

#RSCRIPT = '/home/wenlin/local/R-2.15.3/bin/Rscript'
#DELOX = '/home/xcao/p/delox/delox.R'

import os


def run_delox(read1, read2, prefix='delox', adapter = '/home/xcao/p/delox/loxp-adapter.fasta', RSCRIPT = '/home/wenlin/local/R-2.15.3/bin/Rscript', DELOX = '/home/xcao/p/delox/delox.R', rename=True):
    '''
    read1, read2 is the filename of mate pair sequence result
    prefix is the output prefix for delox
    adapter is the fasta file with the sequence for adapters in mate pairs
    run delox rscript
    if rename is True, combine the outfile to 5 output, that is merge all negative and unpaired to prefix+'.single.fastq' 
    '''
    commandline = '{RSCRIPT} {DELOX} {adapter} {read1} {read2} {prefix}'.format(RSCRIPT=RSCRIPT, DELOX=DELOX, adapter=adapter, read1=read1, read2=read2, prefix=prefix)
    os.system(commandline)
    if rename:
        commandline = 'cat {prefix}_read*.negative.fastq {prefix}_read*.unpaired.fastq > {prefix}_read.single.fastq && rm {prefix}_read*.negative.fastq {prefix}_read*.unpaired.fastq'.format(prefix=prefix)
        os.system(commandline)
    print('finish delox')

description = '''
run_delox(read1, read2, prefix='delox', adapter = '/home/xcao/p/delox/loxp-adapter.fasta', rename=True)
    read1, read2 is the filename of mate pair sequence result
    prefix is the output prefix for delox
    adapter is the fasta file with the sequence for adapters in mate pairs
    run delox rscript
    if rename is True, combine the outfile to 5 output, that is merge all negative and unpaired to prefix+'.single.fastq' 
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f','--read1', help = 'read1 of matepair read', required=True)
    parser.add_argument('-r','--read2',help = 'read2 of matepair read', required = True)
    parser.add_argument('-p','--prefix',help = '''output prefix. default="delox" ''',default='delox')
    parser.add_argument('-a','--adapter', help = 'the fasta file with the sequence for adapters in mate pairs, default /home/xcao/p/delox/loxp-adapter.fasta', default = '/home/xcao/p/delox/loxp-adapter.fasta')
    parser.add_argument('-R','--RSCRIPT', help = '''where is the Rscript, default="/home/wenlin/local/R-2.15.3/bin/Rscript"''', default = '/home/wenlin/local/R-2.15.3/bin/Rscript')
    parser.add_argument('-D','--DELOX', help = 'where is DELOX program, default = "/home/xcao/p/delox/delox.R"', default = '/home/xcao/p/delox/delox.R')
    parser.add_argument('-N','--rename', help = 'whether combine and rename single-end reads, default True', default = True)
    f = parser.parse_args()
    run_delox(read1=f.read1, read2=f.read2, prefix=f.prefix, adapter = f.adapter, RSCRIPT = f.RSCRIPT, DELOX = f.DELOX, rename=f.rename)