
from Bio import SeqIO
import os

BPP = '/home2/s185491/p/bpp/bpp-master/bpp'
TREE1 = '(((J_c_coenia,(J_c_grisea,J_nigrosuffusa)),J_nigrosuffusaTX),(J_neildi,J_z_zonalis));'
TREE2 = '((((J_c_coenia,J_c_grisea),J_nigrosuffusa),J_nigrosuffusaTX),(J_neildi,J_z_zonalis));'
TREE3 = '(((J_c_coenia,J_nigrosuffusaTX),(J_c_grisea,J_nigrosuffusa)),(J_neildi,J_z_zonalis));'
TREES = [TREE1,TREE2,TREE3]
#FOLDER_WORKING = '/work/biophysics/s185491/2018junonia/20181105bpp/temp/'
#FOLDER_OUTPUT1 = '/work/biophysics/s185491/2018junonia/20181105bpp/100frag_bppBPP1/'
#FOLDER_OUTPUT2 = '/work/biophysics/s185491/2018junonia/20181105bpp/100frag_bppBPP2/'
#FOLDER_OUTPUT3 = '/work/biophysics/s185491/2018junonia/20181105bpp/100frag_bppBPP3/'
#FOLDER_OUTPUTS = [FOLDER_OUTPUT1,FOLDER_OUTPUT2,FOLDER_OUTPUT3]
IMAPFILE = '/work/biophysics/s185491/2018junonia/20181105bpp/Imap'
THETA_A = 3
THETA_B = 0.008
TAU_A = 3
TAU_B = 0.04

A01 = '''
          seed =  {SEED}

       seqfile = {SEQFILE}
      Imapfile = {IMAPFILE}
       outfile = {OUTFILE}
      mcmcfile = {MCMCFILE}

  speciesdelimitation = 0 * fixed species tree
         speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 6  J_c_coenia  J_c_grisea  J_neildi  J_nigrosuffusa J_nigrosuffusaTX J_z_zonalis
                    1  1  1  1  1  1
                 {TREE}
        diploid =   0  0  0  0  0  0
   
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 1  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = {THETA_A} {THETA_B}   # invgamma(a, b) for theta
      tauprior = {TAU_A} {TAU_B}    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 40000
      sampfreq = 20
       nsample = 1000000

'''

def alignmentFasta2PhylipBPP(filename):
    '''
    filename is a filename of alignments in fasta format
    if outfilename is None, write to filename+'.phylip'
    '''
    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    workingFolder = os.path.dirname(dirname)
    FOLDER_WORKING = os.path.join(workingFolder,'temp/')
    if not os.path.exists(FOLDER_WORKING):
        os.makedirs(FOLDER_WORKING)
    FOLDER_OUTPUTS = [dirname+'BPP%d/'%e for e in [1,2,3]]
    for f in FOLDER_OUTPUTS:
        if not os.path.exists(f):
            os.makedirs(f)
    outfilename = os.path.join(FOLDER_WORKING,basename + '.phylip')
    fout = open(outfilename,'w')
    
    seqNum = 0
    seqIDs = []
    for s in SeqIO.parse(filename,'fasta'):
        seqNum += 1
        seqIDs.append(s.id)
    
    nameLen = max(len(e) for e in seqIDs) + 3
    
    seqLen = len(s.seq)
    fout.write('%d %d\n'%(seqNum,seqLen))
    
    for s in SeqIO.parse(filename,'fasta'):
        fout.write(s.id + '^'+ s.id +' '*(nameLen - len(s.id))+str(s.seq)+'\n')
    fout.close()
    print('Format converting Finished!')
    
    # make A01 file and get commandline
    commandlines = []
    for n in range(3):
        outfile_A01 = os.path.join(FOLDER_WORKING,basename+'.A01_%d'%(n+1))
        outfile_result = os.path.join(FOLDER_WORKING, basename+'.out_%d'%(n+1))
        outfile_mcmc = os.path.join(FOLDER_WORKING, basename+'.mcmc_%d'%(n+1))
        fout = open(outfile_A01,'w')
        fout.write(A01.format(SEED = -1, SEQFILE=outfilename, IMAPFILE=IMAPFILE, OUTFILE = outfile_result, MCMCFILE = outfile_mcmc, TREE = TREES[n],THETA_A = THETA_A, THETA_B = THETA_B, TAU_A = TAU_A, TAU_B = TAU_B))
        fout.close()
        commandline = '{BPP} --cfile {outfile_A01} && cp -f {outfile_result} {FOLDER_OUTPUT}{basename}'.format(BPP=BPP, outfile_A01 = outfile_A01, outfile_result=outfile_result, FOLDER_OUTPUT = FOLDER_OUTPUTS[n], basename=basename)
        commandlines.append(commandline)
    
    commandlines.append('wait')
    os.system(' & '.join(commandlines))
    print('Finish BPP with tree')
    



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a fasta file of sequence alignments, run BPP for 6 merged genome fragments 100kb, each with 3 individual runs.')
    parser.add_argument('-i','--input', help = 'input file fasta file', required=True)
    f = parser.parse_args()
    alignmentFasta2PhylipBPP(filename=f.input)