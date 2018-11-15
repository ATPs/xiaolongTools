
from Bio import SeqIO

def alignmentFasta2Phylip(filename,outfilename=None):
    '''
    filename is a filename of alignments in fasta format
    if outfilename is None, write to filename+'.phylip'
    '''
    if outfilename is None:
        outfilename = filename + '.phylip'
    fout = open(outfilename,'w')
    
    seqNum = 0
    seqIDs = []
    for s in SeqIO.parse(filename,'fasta'):
        seqNum += 1
        seqIDs.append(s.id)
    
    nameLen = max(len(e) for e in seqIDs) + 3
    
    seqLen = len(s.seq)
    fout.write('%d %d\n'%(seqNum,seqLen))
    
    for s in SeqIO.parse(filename,'fasta'):
        fout.write(s.id + ' '*(nameLen - len(s.id))+str(s.seq)+'\n')
    fout.close()
    print('done!')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a fasta file of sequence alignments, output a phylip file.')
    parser.add_argument('-i','--input', help = 'input file fasta file', required=True)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored',default=None)
    f = parser.parse_args()
    alignmentFasta2Phylip(filename=f.input, outfilename=f.output)