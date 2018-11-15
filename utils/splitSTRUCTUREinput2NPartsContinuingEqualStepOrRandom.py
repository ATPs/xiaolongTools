
import os
import math

def splitStr2NpartsContinue(sequence,N):
    '''
    sequence is a str, return a list of N substr, with almost the same size
    the parts and each parts with continuing bases. 
    eg: 12345678, N=2, return 1234 5678
    '''
    seqlen = len(sequence)
    step = math.ceil(seqlen / N)
    sequences = [sequence[i:i+step] for i in range(0,seqlen,step)]
    return sequences

def splitStr2NpartsEqualStep(sequence, N):
    '''
    sequence is a str, return a list of N substr, with almost the same size
    the parts each parts with bases with the same gap. 
    eg: 12345678, N=2, return 1357 2468
    '''
    sequences = []
    for i in range(N):
        sequences.append([])
    for i,base in enumerate(sequence):
        sequences[i % N].append(base)
    return sequences


def splitStr2NpartsRandomChoice(sequence, N, seed=0):
    '''
    sequence is a str, return a list of N substr, with almost the same size
    the parts and continuing parts. 
    eg: 12345678, N=2, return 1348 7256
    '''
    import random
    random.seed(seed)
    random.shuffle(sequence)
    sequences = splitStr2NpartsContinue(sequence,N)
    return sequences

def splitAlignments2Nparts(filename, N=100, start=1, outfolder='.', seed = 0, method = 0):
    '''
    filename is a alignment in fasta format. Split it to N alignments, each parts with almost the same number of bases
    output the files to outfolder. outfilename will be filename plus '.splitN.0/1/2/...'
    '''
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    fastaname = os.path.basename(filename)
    seqs = open(filename,'r')
    outfilenames = [os.path.join(outfolder,fastaname+'.split'+str(N)+'.'+str(i)) for i in range(N)]
    for outfilename in outfilenames:
        if os.path.exists(outfilename):
            os.remove(outfilename)
            print(outfilename,'removed to store new file')
    outfiles = [open(e,'a') for e in outfilenames]
    
    for seq in seqs:
        seq_elements = seq.split()
        sequence = seq_elements[start:]
        seq_id = seq_elements[:start]
        if method == 0:
            sequences = splitStr2NpartsContinue(sequence,N)
        elif method == 1:
            sequences = splitStr2NpartsEqualStep(sequence,N)
        elif method == 2:
            sequences = splitStr2NpartsRandomChoice(sequence,N,seed)
        else:
            print('wrong method!')
            quit()
        for i in range(len(sequences)):
            outfiles[i].write('\t'.join(seq_id+sequences[i])+'\n')
    for f in outfiles:
        f.close()
    print("done!")



description = '''
given a input file for STRUCTURE, split it to N equal parts and output to N files.
split method: 0 for Continue, 1 for EqualStep, 2 for RandomChoice
                eg: 12345678, N=2, return 
                    1234 5678 for Continue
                    1357 2468 for EqualStep
                    3148 7256 for RandomChoice
-S, --start, start split from which column. Default from column 1, and column 0 is sample name.
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input fasta file', required=True)
    parser.add_argument('-o','--output',help = 'folder of the output folder, default ., current folder', required = False, default = '.')
    parser.add_argument('-N','--number', help = 'number of output files. split to N parts', default = 100, type=int)
    parser.add_argument('-S','--start', help = 'start split from which column. columns are 0, 1, ..., default start from column 1',type=int, default = 1)
    parser.add_argument('-s','--seed', help = 'seed for RandomChoice', default = 0, type=int)
    parser.add_argument('-m','--method', help = 'method to split the sequences. three choices, 0 for Continue, 1 for EqualStep, 2 for RandomChoice. default 0', default = 0, type=int, choices = [0,1,2])
    f = parser.parse_args()
    splitAlignments2Nparts(filename=f.input, N=f.number, outfolder=f.output, seed=f.seed, method=f.method, start=f.start)