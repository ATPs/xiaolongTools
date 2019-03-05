folder = '/home/xcao/w/20180905Junonia_coenia/20181004mitochondrion/'
outfile = '/home/xcao/w/20180905Junonia_coenia/20190131New/mito/20190219Junonia_mito.fa'

from Bio import SeqIO
import os

def collectFa(filename):
    '''
    filename is a assemble result of mitochondria
    return a txt for writing, in fasta format
    '''
    seqkeep = []
    total_len = 0
    for s in SeqIO.parse(filename,'fasta'):
        if len(s.seq) > 500:
            seqkeep.append(s)
            total_len += len(s.seq)
    
    if len(seqkeep) == 0 or total_len < 10000:
        return ''
    
    sample_name = os.path.basename(filename).split('_')[0]
    ls_txt = []
    for n,s in enumerate(seqkeep):
        ls_txt.append('>' + sample_name + '_' + str(n) + '\n' + str(s.seq) +'\n')
    
    return ''.join(ls_txt)

files = os.listdir(folder)
fout = open(outfile,'w')
for f in files:
    fout.write(collectFa(os.path.join(folder,f)))
fout.close()

#get mitochondria sequence based on reference mito
filename = '/home/xcao/w/20180905Junonia_coenia/20190131New/mito/20190219Junonia_mito.fa'
reference = '/home/xcao/w/20180905Junonia_coenia/20190131New/mito/reference.fa'

from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import io
import copy
from collections import Counter
import re

def muscleAlignment(seqs, muscle_exe = "muscle"):
    '''
    align sequences with muscle
    given a list of seqs in SeqIO format,
    return a aligned seqs in SeqIO format
    '''
    f_mem = io.StringIO()
    SeqIO.write(seqs,f_mem,'fasta')
    data = f_mem.getvalue()
    muscle_cline = MuscleCommandline(muscle_exe)
    stdout, stderr = muscle_cline(stdin=data)
    return list(SeqIO.parse(io.StringIO(stdout),'fasta'))

def getAlignLen(seq, minFrag = 10):
    '''
    seq is like 'ATAAA---A-GG-ATCG'
    return alignment length
    do not count fragments shorter than minFrag
    '''
    l = re.split('-+',seq)
    l = [e for e in l if len(e) >= minFrag]
    return sum([len(e) for e in l])

def getAlign(seq1,refseq):
    '''
    seq1 is a Bio.Seq object, refseq is a reference Bio.Seq object
    return str with the same length of refseq, based on alignment result of muslce
    refseq can form a circle, so need to make a new two concatenated copied of reference
    '''
    seqlen = len(refseq.seq)
    refstr = str(refseq.seq).upper()
    r = copy.deepcopy(refseq)
    r.seq = r.seq+r.seq
    aligns = muscleAlignment([seq1, r])
    s1 = str(aligns[0].seq).upper()
    s2 = str(aligns[1].seq).upper()
    
    seq2 = seq1.reverse_complement()
    aligns = muscleAlignment([seq2, r])
    v1 = str(aligns[0].seq).upper()
    v2 = str(aligns[1].seq).upper()
    
    if getAlignLen(s1) < getAlignLen(v1):
        s1 = v1
        s2 = v2
    
    result = ['-' for e in range(seqlen)]
    N = 0
    for b1, b2 in zip(s1,s2):
        if b2 == '-':
            continue
        if N >= seqlen:
            N = N-seqlen
        if b1 != '-':
            if result[N] == '-':
                result[N] = b1
            elif result[N] != refstr[N]:
                result[N] = b1
        N += 1
    return ''.join(result)
    
    
def getMitoBasedRef(dc_seqs, sample_id, refseq):
    '''
    seqs is a list of SeqIO sequences of one sample
    return a sequence based on refseq, a str in fasta format that is ready to write
    refseq is a mito sequences together
    '''
    seqs = dc_seqs[sample_id]
    refstr = str(refseq.seq).upper()
    seqlen = len(refseq.seq)
    ls_align = [getAlign(e,refseq) for e in seqs]
    result = ['-' for e in range(seqlen)]
    for N, bases in enumerate(zip(*ls_align)):
        bases = [e for e in bases if e != '-']
        bases = Counter(bases).most_common()
        if len(bases) != 0:
            if refstr[N] in [e[0] for e in bases]:
                result[N] = refstr[N]
            else:
                result[N] = bases[0][0]
    seq = ''.join(result)
    print(sample_id,'finished')
    return '>'+sample_id+'\n' + seq+'\n'

seqs = list(SeqIO.parse(filename,'fasta'))
refseq = SeqIO.read(reference,'fasta')

dc_seqs = {}
for s in seqs:
    key = s.id.split('_')[0]
    if key not in dc_seqs:
        dc_seqs[key] = []
    dc_seqs[key].append(s)
print('total samples to analyze',len(dc_seqs))


from multiprocessing import Pool
pool = Pool(48)
results = pool.starmap(getMitoBasedRef,[[dc_seqs, k, refseq] for k in dc_seqs])
pool.close()

outfile = '/home/xcao/w/20180905Junonia_coenia/20190131New/mito/20190219Junonia_mito.align.fa'
open(outfile,'w').write(''.join(results))