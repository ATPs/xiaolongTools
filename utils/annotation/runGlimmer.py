#GLIMMER = '/home/xcao/p/GlimmerHMM/GlimmerHMM/bin/glimmerhmm_linux_x86_64'
#genome = '/home/xcao/w/20181205Hermeuptychia/20181224QianPipleline_3318/18_annotate_denovo/genome.fa'
#traindir = '/home/xcao/w/20181205Hermeuptychia/20181224QianPipleline_3318/17_train_denovo/train_glimmer/3318/'
#output = '/home/xcao/w/20181205Hermeuptychia/20181224QianPipleline_3318/18_annotate_denovo/glimmer.gtf'
import uuid
from Bio import SeqIO
import os
from multiprocessing import Pool
import subprocess


def seq2gff(seq, GLIMMER, traindir):
    '''
    seq is a SeqIO sequence, return txt of it running glimmer
    '''
    tempfile = '/dev/shm/' + str(uuid.uuid4())
    open(tempfile,'w').write('>'+seq.id+'\n'+str(seq.seq)+'\n')
    commandline = '{GLIMMER} {genomeseq} -d {traindir} -n 1 -g'.format(GLIMMER=GLIMMER, genomeseq=tempfile, traindir=traindir)
    result = subprocess.check_output(commandline, shell=True)
    os.remove(tempfile)
    result = result.decode('utf-8')
    return result

def runGlimmerHMM(genome='genome.fa',traindir='traindir',GLIMMER='glimmerhmm', outfile='glimmer.gtf', threads=32):
    '''
    given a genome, run glimmer and output the gff result to outfile
    genome can also be a list of SeqIO elements
    traindir is where the train data is stored
    GLIMMER is the location of GLIMMER program
    outfile is where to store the result
    threads is number of CPUs to run the program
    '''
    if isinstance(genome,str):
        ls_seqs = list(SeqIO.parse(genome,'fasta'))
    elif isinstance(genome, list):
        ls_seqs = genome
    else:
        print('wrong input format for genome')
        return None
    
    fout = open(outfile,'w')
    pool = Pool(threads)
    for n in range(0,len(ls_seqs),threads):
        ls_seqs_torun = ls_seqs[n:n+threads]
        results = pool.starmap(seq2gff, [[e, GLIMMER, traindir] for e in ls_seqs_torun])
        fout.write(''.join(results))
        print('{:.2%} finished, up to {}'.format((n+threads)/len(ls_seqs), n+threads))
    pool.close()
    
    fout.close()
    print('done')

description = '''runGlimmerHMM(genome='genome.fa',traindir='traindir',GLIMMER='glimmerhmm', outfile='glimmer.gtf', threads=32)
    given a genome, run glimmer and output the gff result to outfile
    genome can also be a list of SeqIO elements
    traindir is where the train data is stored
    GLIMMER is the location of GLIMMER program
    outfile is where to store the result
    threads is number of CPUs to run the program
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'genome input, default = "genome.fa"', required=True,default = 'genome.fa')
    parser.add_argument('-d','--traindir',help = 'traindir is the where the train data is stored, default="traindir"', default = "traindir")
    parser.add_argument('-G','--GLIMMER',help = 'where GLIMMER program is stored. default "glimmerhmm"',default="glimmerhmm")
    parser.add_argument('-o','--outfile',help = 'where to store the output gff file. default "glimmer.gtf"',default="glimer.gtf")
    parser.add_argument('-t','--threads',help = 'number of CPUs to run the program. default 32',default=32, type=int)
    f = parser.parse_args()
    runGlimmerHMM(genome=f.genome,traindir=f.traindir,GLIMMER=f.GLIMMER, outfile=f.outfile, threads=f.threads)