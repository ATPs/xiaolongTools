# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 17:23:58 2019

@author: ATPs
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 14:12:11 2018

@author: ATPs
"""
import os
import sys
try:
    sys.path.append('/home2/s185491/p/xiaolongTools/utils')
except:
    sys.path.append('/home/xcao/p/xiaolongTools/utils')
import fasta2nexus
import MrBayes2Newick



TXT_BAYES='''
begin mrbayes;
{startvals}
lset nst={nst} rates={rates};
mcmc Nruns={Nruns} nchains={nchains} ngen={ngen} stoprule=yes stopval={stopval};
sump;
sumt;
end;
'''

TXT_BEGINTREE='''
begin trees;
tree mystarttree = {start_tree}
end;
'''

def runMrBayes(file_in, folder_bayes = None, folder_newick = None, folder_temp = None, MB='/home2/s185491/p/MrBayes/mb', Nruns=2, nchains=16, ngen=1000000, stopval=0.01, start_tree=None, threads=32, nst=1, rates='Equal'):
    '''
    file_in is a input file in fasta format
    convert the file to nexus format
    run MrBayes
    convert the result to newick format
    if folder_bayes = None, folder_newick = None, folder_temp = None
    folder_bayes = folder of file_in +'MrBayes'
    folder_newick = folder of file_in + 'Newick'
    folder_temp = folder file_in + 'Temp'
    MB, location of the execuable mrbayes location
    Nruns, number of runs for mrbayes
    nchains, number of chains for each run
    ngen, number of generation to run maximum
    stopval, when to stop the running
    start_tree, start mrbayes with a starting tree. default None, with no starting tree. it should be a file location
    threads is the mpi threads, default 32
    nst, rates is the setting of lset. Default setting is nst=1, rates=Equal. nst=6, rates=Invgamma is the complex model
    '''
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
    
    from Bio import SeqIO
    seqids = [s.id for s in SeqIO.parse(file_in,'fasta')]
    seqids = set(seqids)
    
    #prepare TXT_BEGINTREE
    if start_tree is None:
        txt_begintree = ''
        startvals = ''
    else:
        start_tree = open(start_tree).read()
        start_tree = start_tree.split()
        if len(start_tree) == 1:
            start_tree= start_tree[0]
        elif len(start_tree) == 2:
            start_tree = start_tree[1]
        if not start_tree.endswith(';'):
            start_tree = start_tree+';'
        from ete3 import Tree
        t = Tree(start_tree)
        for l in t.iter_leaves():
            if l.name not in seqids:
                l.delete()
        start_tree = t.write(format=9)
        txt_begintree = TXT_BEGINTREE.format(start_tree=start_tree)
        startvals = "startvals tau = mystarttree;"
        open(file_nex,'a').write(txt_begintree)
    
    #prepare TXT_BAYES
    txt_bayes = TXT_BAYES.format(Nruns=Nruns, nchains=nchains, ngen=ngen, stopval=stopval,startvals=startvals, rates=rates, nst=nst)
    # append MrBayes commands
    open(file_nex,'a').write(txt_bayes)
    
    # run MrBayes
    os.system('module load openmpi/intel')
    os.system('module load beagle-lib')
    
    if threads > 1:
        commandline = 'cd {folder_temp} && mpirun -np {threads} --oversubscribe {MB} {file_nex}'.format(MB=MB, file_nex=basename,threads=threads, folder_temp=folder_temp)
    elif threads == 1:
        commandline = 'cd {folder_temp} && {MB} {file_nex}'.format(MB=MB, file_nex=basename,threads=threads, folder_temp=folder_temp)
    else:
        print('wrong threads')
        return None
    os.system(commandline)
    
    # put result to folder_bayes
    file_bayes = os.path.join(folder_bayes, basename)
    os.rename(file_nex+'.con.tre',file_bayes)
    
    # convert output of MrBayes to newick format
    file_newick = os.path.join(folder_newick,basename)
    MrBayes2Newick.readMrBayesTree(file_in=file_bayes, file_out=file_newick)
    print('Done',basename)
    return None

description = '''
    file_in is a input file in fasta format
    convert the file to nexus format
    run MrBayes
    convert the result to newick format
    if folder_bayes = None, folder_newick = None, folder_temp = None
    folder_bayes = folder of file_in +'MrBayes'
    folder_newick = folder of file_in + 'Newick'
    folder_temp = folder file_in + 'Temp'
    MB, location of the execuable mrbayes location
    Nruns, number of runs for mrbayes
    nchains, number of chains for each run
    ngen, number of generation to run maximum
    stopval, when to stop the running
    start_tree, start mrbayes with a starting tree. default None, with no starting tree. it should be a file location
    threads is the mpi threads, default 32
    nst, rates is the setting of lset. Default setting is nst=1, rates=Equal. nst=6, rates=Invgamma is the complex model
'''
if __name__ == '__main__':
    print(description)
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--file_in', help = 'input file of fasta file', required=True)
    parser.add_argument('-b','--folder_bayes', help = 'where to store the result mrbayes nexus tree, default None, file_in folder.MrBayes', default=None)
    parser.add_argument('-n','--folder_newick', help = 'where to store the result newick file, default None, file_in folder.Newick', default=None)
    parser.add_argument('-T','--folder_temp', help = 'temp folder to running mrbayes, default None, file_in folder.Temp', default=None)
    parser.add_argument('-M','--MB', help = 'location of mrbayes, default for bioHPC:/home2/s185491/p/MrBayes/mb', default = '/home2/s185491/p/MrBayes/mb')
    parser.add_argument('-N','--Nruns', help = 'number of indenpendent runs for mrbayes, default 2', default=2, type=int)
    parser.add_argument('-c','--nchains', help = 'number of chains for each run in mrbayes, default 16', default=16, type=int)
    parser.add_argument('-g','--ngen', help = 'number of generations maximum for mrbayes, default 1000000', default=1000000, type=int)
    parser.add_argument('-s','--stopval', help = 'stop value for mrbayes run, default 0.01', default=0.01, type=float)
    parser.add_argument('-S','--start_tree', help = 'location of starting tree in newick format, default None', default=None)
    parser.add_argument('-t','--threads', help = 'number of CPU for mpirun, default 32', default=32, type=int)
    f = parser.parse_args()
    runMrBayes(file_in=f.file_in, folder_bayes = f.folder_bayes, folder_newick = f.folder_newick, folder_temp = f.folder_temp, MB=f.MB, Nruns=f.Nruns, nchains=f.nchains, ngen=f.ngen, stopval=f.stopval, start_tree=f.start_tree, threads=f.threads)