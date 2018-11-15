
from Bio import SeqIO
from multiprocessing import Pool


def calculateT1T2_4sequencesOrdered(seq1, seq2, seq3, seq4):
    '''
    the phylogenetc topology is "((seq1,seq2),seq3),seq4", where seq4 is a outgroup
    seq4 defines the ancestral allele (A), 
    T2 = (n_ABAA + n_BAAA)/2N
    T1 = T2 + n_BBAA/N
    where N is the total number of sites in the sampled region.
    For this function, gap in the sequence is not allowed. Remove gaps before use
    return T1 and T2
    '''
    N = len(seq1)
    n_ABAA = 0
    n_BAAA = 0
    n_BBAA = 0
    for P1,P2,P3,P4 in zip(seq1,seq2,seq3,seq4):
        if P1 == P4 and P2 != P4 and P3 == P4:
            n_ABAA += 1
        elif P1 != P4 and P2 == P4 and P3 == P4:
            n_BAAA += 1
        elif P1 != P4 and P2 != P4 and P3 == P4:
            n_BBAA += 1
    T2 = (n_ABAA + n_BAAA) /(2*N)
    T1 = T2 + n_BBAA/N
    return T1, T2

def removeGapsListSeq(seqs):
    '''
    seqs is a list of sequences with equal length.
    remove sites with gap symbols. 
    return a list of the filtered sequences.
    '''
    bases = 'ATCG'
    seqs_keep = []
    for sites in zip(*seqs):
        if all([e in bases for e in sites]):
            seqs_keep.append(sites)
    seqs_final = []
    for seq in zip(*seqs_keep):
        seqs_final.append(''.join(seq))
    
    return seqs_final

def calculateT1T2FromFastaAlignment(seqs_fasta, orders=None,minlen=0):
    '''
    seqs_fasta is a filename with the fasta sequences, or a list of fasta sequences in SeqIO format
    orders is a filename with seq1, seq2, seq3, seq4 in correct order (split by space or tab), or a list with ids of them. If orders is None, then we use the original order of sequences in seqs_fasta as the input for calculateT1T2_4sequencesOrdered
    if the filtered seqs is shorter than minlen, return None
    '''
    if type(seqs_fasta) is str:
        seqs_fasta = list(SeqIO.parse(seqs_fasta,'fasta'))
    
    if orders is None:
        seqs = [str(s.seq) for s in seqs_fasta]
        if len(seqs) != 4:
            print('warning, calculateT1T2FromFastaAlignment best works with 4 sequences. total sequences are', len(seqs))
    else:
        dc_seqs = {s.id:str(s.seq) for s in seqs_fasta}
        if type(orders) is str:
            orders = open(orders).read().split()
        if type(orders) is not list:
            print('orders not None, not list, something wrong')
            return None
        if len(orders) != 4:
            print('calculateT1T2FromFastaAlignment only works with 4 sequences.')
            return None
        seqs = [dc_seqs[e] for e in orders]
    
    seqs_filtered = removeGapsListSeq(seqs)
    
    seqlen = len(seqs_filtered[0])
    if seqlen < minlen:
        print('after filtering, sequences have length of {seqlen}, which is shorter than {minlen}'.format(seqlen=seqlen,minlen=minlen))
        return None
    
    return calculateT1T2_4sequencesOrdered(*seqs_filtered)

def calculateT1T2FromMultipleFastaAlignment(file_fastas, ls_orders, minlen = 0, threads = 16, output=None):
    '''
    file_fastas is a file storing the location of aligned fasta files, or a list with the locations
    similar to orders in calculateT1T2FromFastaAlignment, ls_orders is a list of orders that work with each file in file_fastas. In each element, the last one is the outgroup. ls_orders can be a list of list, or a filename, with each line for one fasta file
    minlen is the minimum length of sequences after removing sites with gaps to calculate T1 and T2
    threads is the threads to run in parrallel
    if output is not None, write the result to output, in the format of
        filename_of_fasta: T1, T2
    return a list with filename_of_fasta as key and (T1,T2) as value
    '''
    if type(file_fastas) is str:
        file_fastas = open(file_fastas).readlines()
        file_fastas =[e.strip() for e in file_fastas]
    elif type(file_fastas) is not list:
        print('file_fastas not a filename, not a list. Wrong!')
        return None
    
    if type(ls_orders) is str:
        ls_orders = open(ls_orders).readlines()
        ls_orders = [e.split() for e in ls_orders]
    elif type(ls_orders) is not list:
        print('''target_samples not a filename, not in 'sample1,sample2' format, not a list. Wrong!''')
        return None
    
    pool = Pool(threads)
    results = pool.starmap(calculateT1T2FromFastaAlignment,[[file_fastas[e], ls_orders[e],minlen] for e in range(len(file_fastas))])
    pool.close()
    
    if output is not None:
        with open(output,'w') as f:
            for k,v in zip(file_fastas,results):
                if v is not None:
                    f.write('{k}:{v0},{v1}\n'.format(k=k,v0=v[0],v1=v[1]))
    
    dc_results = {k:v for k,v in zip(file_fastas,results)}
    return dc_results

description = '''file_fastas is a file storing the location of aligned fasta files, or a list with the locations
    similar to orders in calculateT1T2FromFastaAlignment, ls_orders is a list of orders that work with each file in file_fastas. In each element, the last one is the outgroup. ls_orders can be a list of list, or a filename, with each line for one fasta file
    the phylogenetc topology is "((seq1,seq2),seq3),seq4", where seq4 is a outgroup
    seq4 defines the ancestral allele (A), 
    T2 = (n_ABAA + n_BAAA)/2N
    T1 = T2 + n_BBAA/N
    minlen is the minimum length of sequences after removing sites with gaps to calculate T1 and T2
    threads is the threads to run in parrallel
    if output is not None, write the result to output, in the format of
        filename_of_fasta: T1, T2
    return a list with filename_of_fasta as key and (T1,T2) as value
'''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'file_fastas is a file storing the location of aligned fasta files, or a list with the locations', required=True)
    parser.add_argument('-t','--threads',help = 'number of threads to use, default = 16', default = 16,type=int)
    parser.add_argument('-l','--orders',help = 'list of ordered four samples for each input fasta file. the phylogenetc topology is "((seq1,seq2),seq3),seq4"', required = True)
    parser.add_argument('-o','--output',help = 'where to store the output. default=None, only return the result', default = None)
    parser.add_argument('-m','--minlen', help = '', default = 0, type=int)
    f = parser.parse_args()
    calculateT1T2FromMultipleFastaAlignment(file_fastas = f.input, ls_orders = f.orders, minlen = f.minlen, threads = f.threads, output=f.output)
    
