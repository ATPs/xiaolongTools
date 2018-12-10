# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 14:12:11 2018

@author: ATPs
"""
import os
import sys
sys.path.append('/home2/s185491/p/xiaolongTools/utils')
import fasta2nexus
import MrBayes2Newick



TXT_BAYES='''begin mrbayes;
mcmc Nruns=2 nchains=16 ngen=2000 samplefreq=10 stoprule=yes stopval=0.05  Nswaps=3;
sump;
sumt;
end;'''
MB = '/home2/s185491/p/MrBayes/mb'

def runMrBayes(file_in, folder_bayes = None, folder_newick = None, folder_temp = None):
    '''
    file_in is a input file in fasta format
    convert the file to nexus format
    run MrBayes
    convert the result to newick format
    if folder_bayes = None, folder_newick = None, folder_temp = None
    folder_bayes = folder of file_in +'MrBayes'
    folder_newick = folder of file_in + 'Newick'
    folder_temp = folder file_in + 'Temp'
    '''
    os.system('module load openmpi/intel/3.1.1')
    os.system('module load beagle-lib/2.1.2')
    # prepare folders
    folder_in = os.path.dirname(file_in)
    basename = os.path.basename(file_in)
    if folder_bayes is None:
        folder_bayes = folder_in + 'MrBayes'
    if folder_newick is None:
        folder_newick = folder_in + 'Newick'
    if folder_temp is None:
        folder_temp = folder_in + 'Temp'
    if not os.path.exists(folder_bayes):
        os.makedirs(folder_bayes)
    if not os.path.exists(folder_newick):
        os.makedirs(folder_newick)
    if not os.path.exists(folder_temp):
        os.makedirs(folder_temp)
    
    # convert fasta to nexus format
    file_nex = os.path.join(folder_temp,basename)
    fasta2nexus.changeAlignment2Nex(file_in=file_in, file_out=file_nex)
    
    # append MrBayes commands
    open(file_nex,'a').write(TXT_BAYES)
    
    # run MrBayes
    commandline = 'mpirun -np 32 --oversubscribe {MB} {file_nex}'.format(MB=MB, file_nex=file_nex)
    os.system(commandline)
    
    # put result to folder_bayes
    file_bayes = os.path.join(folder_bayes, basename)
    os.rename(file_nex+'.con.tre',file_bayes)
    
    # convert output of MrBayes to newick format
    file_newick = os.path.join(folder_newick,basename)
    MrBayes2Newick.readMrBayesTree(file_in=file_bayes, file_out=file_newick)
    print('Done',basename)
    return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a fasta file, run MrBayes and output MrBayes result and newick result')
    parser.add_argument('-i','--input', help = 'input file of fasta file', required=True)
    parser.add_argument('-b','--folder_bayes', help = 'input file of fasta file', default=None)
    parser.add_argument('-n','--folder_newick', help = 'input file of fasta file', default=None)
    parser.add_argument('-t','--folder_temp', help = 'input file of fasta file', default=None)
    f = parser.parse_args()
    runMrBayes(file_in=f.input,folder_bayes = f.folder_bayes, folder_newick = f.folder_newick, folder_temp = f.folder_temp)