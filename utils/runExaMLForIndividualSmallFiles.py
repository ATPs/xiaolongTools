import os
import shutil

ExaML_parser = '/home2/s185491/p/examl/ExaML/parser/parse-examl'
ExaML = '/home2/s185491/p/examl/ExaML/examl/examl'
mpirun = '/cm/shared/apps/openmpi/intel/3.1.1/bin/mpirun'

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


def runExaMLForIndividualSmallFiles(filename, refTree, threads=1):
    '''
    filename is the full path of a aligned fasta file
    output and save the best tree in a folder by adding '.ExaMLbestTree' to the folder name.
    filename = '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50/exon78368'
    bestTree is stored in '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50.ExaMLbestTree/exon78368'
    '''
    folder = os.path.dirname(filename)
    name = os.path.basename(filename)
    outfolder = folder + '.ExaMLbestTree/'
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    
    basefolder = os.path.basename(folder)
    workFolder = '/dev/shm/'+basefolder + name
    if not os.path.exists(workFolder):
        os.makedirs(workFolder)
    os.chdir(workFolder)
    
    convertFasta2ExaMLinput(filename,outfile=name)
    
    commandline = '{ExaML_parser} -s {name} -n {name} -m DNA  '.format(name=name,ExaML_parser=ExaML_parser)
    os.system(commandline)
    os.system('module add openmpi/intel/3.1.1')
    commandline = '{mpirun} -np {threads} {ExaML} -m GAMMA -t {refTree} -s {name}.binary -n {name}'.format(mpirun=mpirun, ExaML=ExaML, threads=threads, name=name, refTree = refTree)
    os.system(commandline)
    shutil.copy('ExaML_result.' + name, outfolder + name)
    shutil.rmtree(workFolder)
    print('done for', name)

description = '''
    filename is the full path of a aligned fasta file
    output and save the best tree in a folder by adding '.ExaMLbestTree' to the folder name.
    filename = '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50/exon78368'
    bestTree is stored in '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50.ExaMLbestTree/exon78368'
    intermediate file stored in memory. Thus, only good for small files.
    default, use all CPUs
    -r, --refTree, the reference tree required by ExaML
    '''
if __name__ == '__main__':
    import argparse
    import multiprocessing
    CPU_COUNTS = multiprocessing.cpu_count()
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file storing the location of aligned fasta file', required=True)
    parser.add_argument('-r','--refTree', help = 'input reference Tree to guid ExaML', required=True)
    parser.add_argument('-T', '--threads',help = 'number of threads to run ExaML', default = CPU_COUNTS)
    f = parser.parse_args()
    runExaMLForIndividualSmallFiles(f.input, refTree=f.refTree, threads=f.threads)
