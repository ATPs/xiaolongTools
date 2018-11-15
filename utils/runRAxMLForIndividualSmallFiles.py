import os
import shutil
import sys

filename = sys.argv[1]

def runRAxMLForIndividualSmallFiles(filename, threads=1):
    '''
    filename is the full path of a aligned fasta file
    output and save the best tree in a folder by adding '.RAxMLbestTree' to the folder name.
    filename = '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50/exon78368'
    bestTree is stored in '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50.RAxMLbestTree/exon78368'
    '''
    folder = os.path.dirname(filename)
    name = os.path.basename(filename)
    outfolder = folder + '.RAxMLbestTree/'
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    
    basefolder = os.path.basename(folder)
    workFolder = '/dev/shm/'+basefolder + name
    os.system('rm -rf /dev/shm/*')
    if os.path.exists(workFolder):
        shutil.rmtree(workFolder)
    if not os.path.exists(workFolder):
        os.makedirs(workFolder)
    os.chdir(workFolder)
    shutil.copy(filename,'./')
    commandline = '/home2/s185491/p/raxml/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 234 -s {name} -n {name} -T {threads}'.format(name=name,threads=threads)
    os.system(commandline)
    shutil.copy('RAxML_bestTree.' + name, outfolder + name)
    shutil.rmtree(workFolder)
    print('done for', name)

description = '''
    filename is the full path of a aligned fasta file
    output and save the best tree in a folder by adding '.RAxMLbestTree' to the folder name.
    filename = '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50/exon78368'
    bestTree is stored in '/work/biophysics/s185491/2018junonia/geneTrees/sitesKeep_individual_exons.rmGAP0.5.minlen50.RAxMLbestTree/exon78368'
    intermediate file stored in memory. Thus, only good for small files
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file storing the location of aligned fasta file', required=True)
    parser.add_argument('-T', '--threads',help = 'number of threads to run RAxML', default = 1)
    f = parser.parse_args()
    runRAxMLForIndividualSmallFiles(f.input, threads=f.threads)
