# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:05:33 2018

@author: ATPs
"""
import os
ExaML_parser = '/home2/s185491/p/examl/ExaML/parser/parse-examl'
ExaML = ['/home2/s185491/p/examl/ExaML/examl/examl','/home2/s185491/p/examl/ExaML/examl/examl-AVX','/home2/s185491/p/examl/ExaML/examl/examl-OMP-AVX']


def convertFasta2ExaMLinput(filename, outfile = None):
    '''
    Given a filename of fasta file, output a file as input of ExaML
    if outfile is None, outfile = filename+'.ExaML'
    '''
    if outfile is None:
        outfile = filename + '.ExaML'
    fout = open(outfile,'w')
    from Bio import SeqIO
    seqs = SeqIO.parse(filename,'fasta')
    s = next(seqs)
    seqlen = len(s.seq)
    n = 1
    for s in seqs:
        n += 1
    print('totally ', n, 'sequences of length', seqlen)
    fout.write(str(n) + ' ' + str(seqlen) +'\n')
    for s in SeqIO.parse(filename,'fasta'):
        fout.write(s.id + ' ' + str(s.seq)+'\n')
    fout.close()

def runExaML(filename, referenceTree=None, threads=320, maxjobsperNode = 0,whichExaml=0):
    '''
    filename is a file of fasta file
    referenceTree is the reference tree to run ExaML
    process the fasta file and run ExaML
    folder is where to store the files. default the same as the input filename
    basename, prefix for output files, default the same as the input filename
    threads, numbers of cpus to use. default 320
    '''
    folder = os.path.dirname(filename)
    basename = os.path.basename(filename)
    
    #convert fasta file to input for ExaML
    outfile = os.path.join(folder,basename+'.ExaML')
    convertFasta2ExaMLinput(filename, outfile)
    
    #generate random tree if referenceTree is None
    if referenceTree is None:
        referenceTree = filename + '.referenceTree'
        from ete3 import Tree
        from Bio import SeqIO
        leafnames = [e.id for e in SeqIO.parse(filename,'fasta')]
        t = Tree()
        t.populate(len(leafnames),names_library=leafnames)
        t.write(outfile=referenceTree)
        
    file_binary = basename + '.binary'
    
    #run ExaML
    os.system('module load openmpi/intel/3.1.1')
    if maxjobsperNode == 0:
        if not os.path.exists(file_binary):
            commandline = 'cd {folder} && {ExaML_parser} -s {basename}.ExaML -n {basename} -m DNA  && mpirun -np {threads} {ExaML} -m GAMMA -t {referenceTree} -s {basename}.binary -n {basename}'.format(ExaML_parser=ExaML_parser,ExaML=ExaML[whichExaml], folder=folder,basename=basename, referenceTree=referenceTree, threads=threads)
        else:
            commandline = 'cd {folder} && mpirun -np {threads} {ExaML} -m GAMMA -t {referenceTree} -s {basename}.binary -n {basename}'.format(ExaML_parser=ExaML_parser,ExaML=ExaML[whichExaml], folder=folder,basename=basename, referenceTree=referenceTree, threads=threads)
    else:
        commandline = 'cd {folder} && {ExaML_parser} -s {basename}.ExaML -n {basename} -m DNA  && mpirun -np {threads} -npernode {maxjobsperNode} {ExaML} -m GAMMA -t {referenceTree} -s {basename}.binary -n {basename}'.format(ExaML_parser=ExaML_parser,ExaML=ExaML[whichExaml], folder=folder,basename=basename, referenceTree=referenceTree, threads=threads,maxjobsperNode=maxjobsperNode)
    print(commandline)
    os.system(commandline)

description = '''
    filename is a file of fasta file
    referenceTree is the reference tree to run ExaML. if not provided, use ete3.Tree to generate a random tree
    process the fasta file and run ExaML
    folder is where to store the files. default the same as the input filename
    basename, prefix for output files, default the same as the input filename
    threads, numbers of cpus to use. default 320
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file storing the location of aligned fasta file', required=True)
    parser.add_argument('-r','--referenceTree', help = 'input reference Tree to guid ExaML. If not provided, will generate a random tree as input', default=None)
    parser.add_argument('-T', '--threads',help = 'number of threads to run ExaML', default = 320, type=int)
    parser.add_argument('-j', '--maxjobsperNode',help = 'maxjobperNode, default 0, no limit', default = 0, type=int)
    parser.add_argument('-e','--whichExaml',help='which examl program to use. 0 for examl, 1 for examl-AVX, 2 for examl-OMP-AVX. default 0',default=0,type=int,choices=[0,1,2])

    f = parser.parse_args()
    runExaML(filename=f.input, referenceTree=f.referenceTree, threads=f.threads, whichExaml=f.whichExaml)