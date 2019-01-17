

import pandas as pd
import glob
from itertools import chain
from multiprocessing import Pool
from Bio import SeqIO


def getLen(ls_qstart_qend):
    '''
    ls_qstart_qend is a list of (qstart, qend)
    return length of unique ids
    '''
    sites = set()
    ls_qstart_qend = set(ls_qstart_qend)
    for qstart, qend in ls_qstart_qend:
        sites = sites.union(set(list(range(int(qstart),int(qend)+1))))
    return len(sites)


def improve_genome(file_genome='genome', file_stats='merged.stats', file_blast = 'all.blast',depth_min=20, outfile='genome.improved', identity_min = 95, threads=16, maxMissing=200):
    '''
    file_genome is the genome to imporve
    file_stats is the file generated from get_coverage_histogram, with coverage of each scaffolds
    depth_min is the minimum required depth to keep a scaffold
    outfile is the where to store the improved genome
    identity_min is the filter for all.blast to compare two scaffolds
    file_blast is the output of blast the genome against itself. accept "*" in file_blast name for multiple files
    blast is run as 
        blastn -query wenlin3303min200.fa -db wenlin3303min200.fa -out all.blast -outfmt '6 qseqid sseqid pident evalue qlen qstart qend slen sstart send' -num_threads 48 -word_size 50
    modified from Qian and Jing's code
    '''
    st_scf_genome = set()
    for s in SeqIO.parse(file_genome,'fasta'):
        st_scf_genome.add(s.id)
    print(len(st_scf_genome),'scaffolds in input genome file')
    
    df = pd.read_csv(file_stats,sep='\t')
    print(df.shape[0],'scaffolds in file_stats')
    df =  df[df['scf'].apply(lambda x:x in st_scf_genome)]
    scf_keep = set(df[df['depth'] >= depth_min]['scf'])
    scf_checking = set(df[df['depth'] < depth_min]['scf'])
    dc_scflen = dict(zip(df['scf'], df['scf_len_noN']))
    print('finish processing stats')
    print(len(scf_keep),'scaffold will be kept')
    print(len(scf_checking), 'scaffolds need to be checked')
    
    dc_scfcover = {}
    file_blasts = glob.glob(file_blast)
    ls_fo = [open(f) for f in file_blasts]
    fo = chain(*ls_fo)
    for line in fo:
        qseqid,sseqid,pident,evalue,qlen,qstart,qend,slen,sstart,send = line.split()
        if qseqid in scf_checking and sseqid in scf_keep and qseqid != sseqid and float(pident) > identity_min:
            if qseqid not in dc_scfcover:
                dc_scfcover[qseqid] = []
            dc_scfcover[qseqid].append((qstart, qend))
    print('finish scanning blast file')
    
    keys = dc_scfcover.keys()
    values = dc_scfcover.values()
    pool = Pool(threads)
    ls_len = pool.map(getLen, list(values))
    pool.close()
    dc_scfcover = dict(zip(keys,ls_len))
    
    for scf in dc_scfcover:
        #dc_scfcover[scf] = len(dc_scfcover[scf])
        qcover = dc_scfcover[scf]
        qlen = dc_scflen[scf]
        cover_rate = qcover / qlen
        uncover_len = max(qlen -qcover, 0) #in case qlen < qcover, set unconver_len to 0
        if cover_rate/0.9 > uncover_len/maxMissing and uncover_len < maxMissing:
            continue
        else:
            scf_keep.add(scf)
    print('finish filtering')
    print(len(scf_keep),'scaffolds kept')
    
    scf_total_count = 0
    scf_total_len = 0
    keep_total_count = 0
    keep_total_len = 0
    fout = open(outfile,'w')
    for s in SeqIO.parse(file_genome,'fasta'):
        seqlen = len(s.seq)
        scf_total_count += 1
        scf_total_len += seqlen
        if s.id in scf_keep:
            keep_total_count += 1
            keep_total_len += seqlen
            fout.write('>'+s.description+'\n'+str(s.seq)+'\n')
    fout.close()
    print('before filtering, there are', scf_total_count,'scaffolds with total length of',scf_total_len)
    print('after filtering, there are', keep_total_count, 'scaffolds with total length of',keep_total_len)
    print('%.4f scaffold left, %.4f bases left'%(keep_total_count/scf_total_count, keep_total_len/scf_total_len))
    print('done, finish writing output\n')


description = '''file_genome is the genome to imporve
    file_stats is the file generated from get_coverage_histogram, with coverage of each scaffolds
    depth_min is the minimum required depth to keep a scaffold
    outfile is the where to store the improved genome
    identity_min is the filter for all.blast to compare two scaffolds
    file_blast is the output of blast the genome against itself. accept "*" in file_blast name for multiple files
    blast is run as 
        blastn -query wenlin3303min200.fa -db wenlin3303min200.fa -out all.blast -outfmt '6 qseqid sseqid pident evalue qlen qstart qend slen sstart send' -num_threads 48 -word_size 50
    modified from Qian and Jing's code
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--file_genome', help = 'file_genome is the genome to imporve, default="genome" ', default='genome')
    parser.add_argument('-s','--file_stats', help = 'file_stats is the file generated from get_coverage_histogram, with coverage of each scaffolds, default="merged.stats" ', default='merged.stats')
    parser.add_argument('-b','--file_blast', help = 'file_blast is the output of blast the genome against itself. default = "all.blast". accept "*" in file_blast name for multiple files', default='all.blast')
    parser.add_argument('-d','--depth_min', help = 'depth_min is the minimum required depth to keep a scaffold, default =20', default=20, type=int)
    parser.add_argument('-o','--outfile', help = 'outfile is the where to store the improved genome, default="genome.improved" ', default='genome.improved')
    parser.add_argument('-m','--identity_min', help = 'identity_min is the filter for all.blast to compare two scaffolds, default=95', default=95, type=int)
    parser.add_argument('-t','--threads',help = 'number of threads to use. default 16', default = 16,type=int)
    parser.add_argument('-M','--maxMissing',help = 'maxMissing bases when remove a scaffold. default 200', default = 200,type=int)
    f = parser.parse_args()
    improve_genome(file_genome=f.file_genome, file_stats=f.file_stats, file_blast = f.file_blast,depth_min=f.depth_min, outfile=f.outfile, identity_min = f.identity_min, threads=f.threads, maxMissing=f.maxMissing)
