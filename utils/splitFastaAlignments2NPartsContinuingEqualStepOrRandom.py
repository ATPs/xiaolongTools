
from Bio import SeqIO
import os

def splitStr2NpartsContinue(sequence,N):
    '''
    sequence is a str, return a list of N substr, with almost the same size
    the parts and each parts with continuing bases. 
    eg: 12345678, N=2, return 1234 5678
    '''
    seqlen = len(sequence)
    step = seqlen // N
    if seqlen % N != 0:
        step = step+1
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
    return [''.join(s) for s in sequences]


def splitStr2NpartsRandomChoice(sequence, N, seed=0):
    '''
    sequence is a str, return a list of N substr, with almost the same size
    the parts and continuing parts. 
    eg: 12345678, N=2, return 1348 7256
    '''
    import random
    random.seed(seed)
    sequence = list(sequence)
    random.shuffle(sequence)
    sequence = ''.join(sequence)
    sequences = splitStr2NpartsContinue(sequence,N)
    return sequences

def splitAlignments2Nparts(filename, N=100, outfolder='.',baseNumberMin=0, seed = 0, method = 0):
    '''
    filename is a alignment in fasta format. Split it to N alignments, each parts with almost the same number of bases
    output the files to outfolder. outfilename will be filename plus '.splitN.0/1/2/...'
    '''
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    fastaname = os.path.basename(filename)
    seqs = SeqIO.parse(filename,'fasta')
    outfilenames = [os.path.join(outfolder,fastaname+'.split'+str(N)+'.'+str(i)) for i in range(N)]
    for outfilename in outfilenames:
        if os.path.exists(outfilename):
            os.remove(outfilename)
            print(outfilename,'removed to store new file')
    outfiles = [open(e,'a') for e in outfilenames]
    
    for seq in seqs:
        sequence = str(seq.seq)
        seq_id = seq.id
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
            if sum(map(lambda x:x in 'ATCG', sequences[i])) > baseNumberMin:
                outfiles[i].write('>'+seq_id+'\n'+sequences[i]+'\n')
    for f in outfiles:
        f.close()
    print("done!")



description = '''
given a alignment in fasta format, split it to N equal parts and output to N files.
some sequences will be removed if the count of "ATCG" bases is less than a baseNumberMin, default 0.
split method: 0 for Continue, 1 for EqualStep, 2 for RandomChoice
                eg: 12345678, N=2, return 
                    1234 5678 for Continue
                    1357 2468 for EqualStep
                    3148 7256 for RandomChoice
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input fasta file', required=True)
    parser.add_argument('-o','--output',help = 'folder of the output folder, default ., current folder', required = False, default = '.')
    parser.add_argument('-N','--number', help = 'number of output files. split to N parts', default = 100, type=int)
    parser.add_argument('-b','--baseNumberMin', help = 'minimum base numbers, default 0', default = 0, type=int)
    parser.add_argument('-s','--seed', help = 'seed for RandomChoice', default = 0, type=int)
    parser.add_argument('-m','--method', help = 'method to split the sequences. three choices, 0 for Continue, 1 for EqualStep, 2 for RandomChoice', default = 0, type=int, choices = [0,1,2])
    f = parser.parse_args()
    splitAlignments2Nparts(filename=f.input, N=f.number, outfolder=f.output, baseNumberMin=f.baseNumberMin, seed=f.seed, method=f.method)