
from Bio import SeqIO
import pandas as pd

def get_coverage_histogram(file_genome='genome', file_readcount='merged.idxstats', readlen=150, outfile_prefix=None):
    '''
    file_genome is the file storing the genome, file_readcount is the file of output of samtools idxstats
    file_readcount
        TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads
    readlen is average read length to calculate coverage
    if outfile_prefix is None, output two file, one is merged.stats, the other is merged.depth_hist
else: outfile_prefix +'.stats'/'.depth_hist'
    '''
#    file_readcount = '/home/xcao/w/20181205Hermeuptychia/20181208QianPipeline/8_improve_assembly/merged.idxstats'
#    file_genome = '/home/xcao/w/genomes/Hermeuptychia/wenlin3303min200.fa'
    
    dc_seqlenNoN = {}
    for s in SeqIO.parse(file_genome,'fasta'):
        seq = str(s.seq)
        seq.upper()
        dc_seqlenNoN[s.id] = len(seq) - seq.count('N')
    
    df = pd.read_csv(file_readcount,header=None,sep='\t')
    df.columns = ['scf','scf_len','mapped_reads','unmapped_reads']
    df = df[df['scf'].apply(lambda x:x in dc_seqlenNoN)]
    df = df.set_index(['scf'])
    df = df[[e != '*' for e in df.index]] #remove the last row, which include summary of unmapped reads
    
    scf_len_noN = [dc_seqlenNoN[e] for e in df.index]
    df['scf_len_noN'] = scf_len_noN
    
    print('sequence with scf_len_noN smaller than 100bp will be removed! they are')
    print(df[df['scf_len_noN']<100])
    df = df[df['scf_len_noN']>=100]
    
    df['depth'] = df['mapped_reads'] * readlen / df['scf_len_noN']
    
    dc_scf2depth = dict(zip(list(df.index),[int(e) for e in list(df['depth'])]))
    dc_depth = {v:0 for v in dc_scf2depth.values()}
    for k,v in dc_scf2depth.items():
        if k in dc_seqlenNoN:
            dc_depth[v] += dc_seqlenNoN[k]
    
#    ls_depth = list(df['depth'])
#    ls_depth = [int(e) for e in ls_depth]
#    from collections import Counter
#    dc_depth = Counter(ls_depth)
    df_depth = pd.DataFrame()
    df_depth['count'] = dc_depth.values()
    df_depth.index = dc_depth.keys()
    df_depth = df_depth.sort_index()
    print(df_depth.iloc[:300])
        
    #save df
    if outfile_prefix is None:
        if file_readcount.endswith('.idxstats'):
            outfile1 = file_readcount.rsplit('.idxstats')[0]+'.stats'
            outfile2 = file_readcount.rsplit('.idxstats')[0]+'.depth_hist'
        else:
            outfile1 = file_readcount +'.stats'
            outfile2 = file_readcount+'.depth_hist'
        df.to_csv(outfile1,sep='\t')
        df_depth.to_csv(outfile2,sep='\t')
    else:
        df.to_csv(outfile_prefix + '.stats',sep='\t')
        df_depth.to_csv(outfile_prefix+'.depth_hist',sep='\t')
    print('done!')

description = '''file_genome is the file storing the genome, file_readcount is the file of output of samtools idxstats
    file_readcount
        TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads
    readlen is average read length to calculate coverage
    if outfile_prefix is None, output two file, one is merged.stats, the other is merged.depth_hist
else: outfile_prefix +'.stats'/'.depth_hist'
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--file_genome', help = 'the file storing the genome. default "genome" ', default='genome')
    parser.add_argument('-r','--file_readcount',help = 'the file of output of samtools idxstats. default "merged.idxstats"', default="merged.idxstats")
    parser.add_argument('-L','--readlen',help = 'the average length of reads', default=150, type=int)
    parser.add_argument('-o','--outfile_prefix',help = '''outfile_prefix is None, output two file, one is merged.stats, the other is merged.depth_hist
else: outfile_prefix +'.stats'/'.depth_hist' ''', default=None)
    f = parser.parse_args()
    get_coverage_histogram(file_genome=f.file_genome, file_readcount=f.file_readcount, readlen=f.readlen, outfile_prefix=f.outfile_prefix)
