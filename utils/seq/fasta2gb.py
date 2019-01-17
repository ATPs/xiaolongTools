from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO


def changeAlignment2Nex(file_in, file_out = None, seq_type='DNA'):
    '''
    convert file_in of fasta format to nexus format
    '''
    if file_out is None:
        file_out = file_in +'.nex'
    if seq_type == 'DNA' or seq_type == 'D':
        alphabet = Gapped(IUPAC.extended_dna)
    elif seq_type == 'protein' or seq_type =='P':
        alphabet = Gapped(IUPAC.extended_protein)
    else:
        print(seq_type, 'is not DNA or protein')
        return None
    
    with open(file_in, "rU") as input_handle:
    
        sequences = list(SeqIO.parse(input_handle, "fasta",alphabet=alphabet))
        fout = open(file_out,'w')
        for s in sequences:
            s.id = s.id.split('|')[0]
            fout.write('>'+s.id+'\n'+str(s.seq)+'\n')
        fout.close()
        sequences = list(SeqIO.parse(file_out, "fasta",alphabet=alphabet))
    
    with open(file_out, "w") as output_handle:
        count = SeqIO.write(sequences, output_handle, "genbank")

    print("Converted %i records" % count)
    return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a fasta file of sequence alignments, output a nexus file.')
    parser.add_argument('-i','--input', help = 'input file fasta file', required=True)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored',default=None)
    parser.add_argument('-t', '--seq_type', help = 'sequence type, DNA or protein. default DNA', choices=['DNA','protein','D','P'], default='DNA')
    f = parser.parse_args()
    changeAlignment2Nex(file_in=f.input, file_out=f.output, seq_type=f.seq_type)