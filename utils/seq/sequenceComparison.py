# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:17:21 2017
try to work with the xml output of blast directly, to avoid using the BLAST program with edited scoring matrix
@author: ATPs
"""



from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import io
from collections import Counter
import time
import os
import pandas as pd
import re
from multiprocessing import Pool



def openfile2lsFasta(filename,fmt = "fasta", removeStar = True):
    """
    given a file name, return a list of fasta in SeqIO format
    """
    from Bio import SeqIO
    l = list(SeqIO.parse(open(filename),fmt))
    if removeStar:
        for s in l:
            if s.seq[-1] == '*':
                s.seq = s.seq[:-1]
    return l

def muscleAlignment(seqs, muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe"):
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

def proteinAlignLength(seqs_aln,mincommon = 10,error_rate = 0.02):
    """
    return alignment length based on the output of proteinPairwiseAlignGlobal
    seq1_aln = "MPKSSSN-DLP"
    seq2_aln = "MPRASSNADLP"
    proteinAlignLength(seq1_aln,seq2_aln,1,0), return 8
    proteinAlignLength(seq1_aln,seq2_aln,3,0), return 6
    proteinAlignLength(seq1_aln,seq2_aln,3,0.5), return 10
    """
    seq1_aln = str(seqs_aln[0].seq)
    seq2_aln = str(seqs_aln[1].seq)
    seqlen = len(seq1_aln)
    if seqlen != len(seq2_aln):
        print("input two seqs are not the same length!")
        return None
    seqcommon = "" #find common elements. if two aa not the same, use #. "MP##SSN-AP"
    for num in range(seqlen):
        if seq1_aln[num] == seq2_aln[num]:
            seqcommon = seqcommon + seq1_aln[num]
        else:
            if seq1_aln[num] == "-" or seq2_aln[num] == "-":
                seqcommon = seqcommon + "-"
            else:
                seqcommon = seqcommon + "#"
    seqs = re.split('-+|#{2,}',seqcommon) #split at gap region or if two amino acids are not the same
    common = 0
    for seq in seqs:
        seqlen = len(seq)
        maxerror = int(seqlen * error_rate)
        error = seq.count("#") 
        if error > maxerror:
            seqlen = seqlen - error +maxerror
        if seqlen < mincommon:
            seqlen = 0
        common += seqlen
#        print(seq,seqlen)
    return common

def getProteinAlignLength(seqs, mincommon = 10, error_rate = 0.02,muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe"):
    '''
    seqs include two sequences in SeqIO format. Return the aligned length by first align sequences 
    with muscleAlignment, then calculate aligned length with proteinAlignLength
    '''
    seqs_aln = muscleAlignment(seqs,muscle_exe = muscle_exe)
    return proteinAlignLength(seqs_aln, mincommon=mincommon,error_rate= error_rate)

def seqs2kmerdic(seqs,kmerlen=20):
    '''
    seqs is a list of seqIO sequences. 
    Return a dictionary, with kmers with of length 20 (kmerlen) as key, 
    value is a list of id of seqs in seqs
    '''
    from collections import defaultdict
    dckmer = defaultdict(set) #kmer dic, kmer with its seqs
    for n, seq in enumerate(seqs):
        seq = str(seq.seq)
        if len(seq) >= kmerlen:
            for i in range(len(seq)+1-kmerlen):
                kmernum = seq[i:i+kmerlen]
                dckmer[kmernum].add(n)
    #change set to list
    for kmernum in dckmer:
        dckmer[kmernum] = list(dckmer[kmernum])
    return dckmer

def getPairsFromTwoListSeqs(seqs1, seqs2, identitymin = 20, max_target = float('inf'),error_rate = 0.02,muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe"):
    '''
    seqs1 and seqs2 are two list of sequences in SeqIO format
    return a list of tuple, seq1_id, seq2_id, 
    for each seq1_id, compare with at most max_target sequences in seqs2
    default, compare all sequences in seqs1 and seqs2 pairs which have at least one identical region with length of 20
    '''
    time1 = time.time()
    df = []
    dc2kmer = seqs2kmerdic(seqs2, kmerlen=identitymin)
    time2 = time.time()
    print('%d seconds, making kmer dictionary for seqs2, kmer length is %d'%(time2 - time1, identitymin))
    for n, seq in enumerate(seqs1):
        s = str(seq.seq)
        #get kmers of s
        skmers = set([s[i:i+identitymin] for i in range(len(s) +1 -identitymin)])
        #get ids in seqs2 which shares a identitymin with s
        s_targets = []
        for kmer in skmers:
            if kmer in dc2kmer:
                s_targets += dc2kmer[kmer]
        if len(s_targets) == 0: # skip if seq have no match in seqs2
            continue
        s_targets = Counter(s_targets)
        s_targets = s_targets.most_common()
        for _n, _target in enumerate(s_targets):
            if _n < max_target:
                df.append((n,_target[0]))
            else:
                break
    time4 = time.time()
    print('%d pairs to compare, total time %d'%(len(df),time4-time1))
    return df

def compare2lsSeqs(seqs1, seqs2, identitymin = 20, outfile = None, max_target = float('inf'),error_rate = 0.02,muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe"):
    '''
    use the output of getPairsFromTwoListSeqs, further calculate matched length
    outfile is a tsv file, with three columns: seqs1_id, seqs2_id, matched_length. There is a headline.
    if outfile exist, first readin all lines of outfile, and write all lines other than the last line back (the last line may be incomplete)
    only calculate matched_length for uncompared pairs.
    return a dataframe with three columns: seqs1_id, seqs2_id, matched_length.
    '''
    df = []
    if isinstance(seqs1,str):
        seqs1 = list(SeqIO.parse(seqs1,'fasta'))
    if isinstance(seqs2,str):
        seqs2 = list(SeqIO.parse(seqs2,'fasta'))
    pairsFinished = []#store the finished pairs
    pairs = getPairsFromTwoListSeqs(seqs1, seqs2, identitymin = identitymin, max_target = max_target,error_rate = error_rate,muscle_exe = muscle_exe)
    print('total pairs to compare:', len(pairs))
    savefile = False
    if outfile is not None:
        savefile = True
        if os.path.isfile(outfile):
            templs = open(outfile).readlines()
            fout = open(outfile,'w')#save back the file
            fout.write(templs[0])
            for _line in templs[1:-1]:
                _line = _line.replace('\n','')
                _es = re.split(',|\t|;| ',_line)
                pairsFinished.append((int(_es[0]),int(_es[1])))
                df.append((int(_es[0]),int(_es[1]),int(_es[2])))
                fout.write('\t'.join(_es)+'\n')
            #check if pairsFinished and pairs the same
            print(len(pairsFinished), 'pairs already finished calculating of matched_length')
        else:
            fout = open(outfile,'w')
            fout.write('seq1_id\tseq2_id\tmatched_length\n')
    pairsFinished = set(pairsFinished)
    if len(pairsFinished - set(pairs))!=0:
        print('the file ', outfile, 'is not from the same setting of this run!')
        return None
    for pair in pairs:
        if pair not in pairsFinished:
            seqList = [seqs1[pair[0]], seqs2[pair[1]]]
            matched_length = getProteinAlignLength(seqList,mincommon=identitymin,error_rate=error_rate,muscle_exe=muscle_exe)
            if savefile:
                fout.write(str(pair[0]) +'\t' +str(pair[1]) +'\t' +str(matched_length)+'\n')
            df.append((pair[0],pair[1],matched_length))
    if savefile:
        fout.close()
    print('df len:', len(df))
    df = pd.DataFrame(df,columns = ['seqs1_id', 'seqs2_id', 'matched_length'])
    return df

def test20170405():
    '''
    test functions above
    '''
    f_fasta1 = r"D:\Insects\ManducaSexta\Interpro\20160612MsSPSPH257MCOTOGS2.txt"
    f_fasta2 = r"D:\Insects\ManducaSexta\Interpro\20170612_MsSPSPH193Annotated.txt"
    seqs1 = openfile2lsFasta(f_fasta1)
    seqs2 = openfile2lsFasta(f_fasta2)
    df = compare2lsSeqs(seqs1,seqs2, outfile = r"D:\P\3Language\Xiaolong\python\list.csv")
    print(df.shape)
    df.loc[:,'seq1name'] = [seqs1[i].id for i in df.loc[:,'seqs1_id']]
    df.loc[:,'seq2name'] = [seqs2[i].id for i in df.loc[:,'seqs2_id']]
    df.loc[:,'seq1len'] = [len(seqs1[i]) for i in df.loc[:,'seqs1_id']]
    df.loc[:,'seq2len'] = [len(seqs2[i]) for i in df.loc[:,'seqs2_id']]
    df.loc[:,'seq1seq'] = [seqs1[i].seq for i in df.loc[:,'seqs1_id']]
    df.loc[:,'seq2seq'] = [seqs2[i].seq for i in df.loc[:,'seqs2_id']]
    df.to_csv('list2.csv')
    open('list.txt','w').write('\n'.join([e.id for e in seqs1]))
    
    
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 07 15:32:37 2014

@author: k
"""


import sys
#sys.path.append("D:\\P\\3Language\\Xiaolong\\python\\XCProject\\fasta\\")

def open_to_list(file_name):
    """
    given the file name including the path
    open a file, and return a list. each list is each line of this file. the newline symbols will be removed
    """
    with open(file_name) as f:
        lines = f.read().splitlines()
    return lines


def count_seq_letters_in_a_row(sequence, letter):
    """
    giving a sequence, count the number of certain letters in a row
    sequence, the sequence you provide
    letter, the letter you want to count in the sequence
    for example, sequence="NNAANNn", letter="N"
    lowercase letter will be changed to uppercase
    return[2,3]
    """
    n=[0]
    m=0
    for i in range(len(sequence)-1):
        if sequence[i].upper()==letter:
            n[m]+=1
            if sequence[i+1].upper()!=letter:
                m+=1
                n.append(0)
    if sequence[-1].upper()==letter:
        n[m]=1
    else:
        n.pop()
    return n

def count_seqs_letters_in_a_row(sequences,letter):
    """
    givin a list of sequence in the imported from Bio.SeqIO.parse(),
    give a list of sequence ID, length, and each letter in a row length
    for example, sequence.seq="NNAANNN", sequence.id="scaffold001",letter="N"
    return letter_num[(sequence.id,length)]=[2,3]
    """
    letter_num={}
    for sequence in sequences:
        sequence_len=len(sequence.seq)
        n=count_seq_letters_in_a_row(sequence.seq, letter)
        if n!=[]:
            letter_num[(sequence.id,sequence_len)]=n
    return letter_num
            
def save_dictionary_with_tuple_list(myfile,mydict):
    """
    save the dictionary generated by count_seqs_letters_in_a_row
    letter_num[(scaffold001,7)]=[2,3]
    saved as "scaffold001\t7\t2\nscaffold001\t7\t3\n"
    """
    for keys in mydict:
        for item_element in mydict[keys]:
            myfile.write(keys[0]+"\t"+str(keys[1])+"\t"+str(item_element)+"\n")

def dic_save(filename, dic_name):
    """
    given a filename, dic_name, and command, if command = "s", save dic to filename
    if command = "r", read filename to dic_name
    """
    fe = open(filename,"w")
    for key in dic_name:
        fe.write(str(key)+"\t"+str(dic_name[key])+"\n")
    fe.close()
    return "Done!"

def dic_read_file_of_list(filename):
    """
    given a filename, read it to dict.
    """
    lines = open_to_list(filename)
    dic_name ={}
    for line in lines:
        dic_name[line.split("\t",1)[0]] = line.split("\t",1)[1]
    return dic_name
            
            



def fasta_length(myfasta,remove_star=True):
    """
    give a list of SeqIO input, return a list based on the length of seq
    return a dictionary of the length of the sequences
    """
    fastalen={}
    for fasta in myfasta:
        if not remove_star:
            fastalen[fasta.id]=len(fasta.seq)
        else:
            if fasta.seq[-1] == "*":
                fastalen[fasta.id]=len(fasta.seq)-1
            else:
                fastalen[fasta.id]=len(fasta.seq)
    return fastalen

def open_fasta_to_list(filename, handle = "fasta"):
    """
    given a file_name, open it with SeqIO, and make a list of fasta seqs
    """
    return list(SeqIO.parse(filename, handle))

def openFasta2NameSeqlist(filename, handle = "fasta", nostar = True):
    """
    open fasta file to two list: name list and seqlist
    namelist is the list of descriptions
    default, remove the star at the end of protein sequence.
    """
    mylist = open_fasta_to_list(filename,handle)
    mynames =[]
    myseqs = []
    for ele in mylist:
        mynames.append(ele.description)
        seq = str(ele.seq)
        if nostar:
            if seq[-1] == "*":
                seq = seq[:-1]
        myseqs.append(seq)
    return mynames, myseqs

def open_selected_fasta_to_list(selected_list, filename, handle = "fasta"):
    """
    given a list to open in a filename, open it with SeqIO, and make a list of fasta seqs
    """
    fasta_list = []
    for sequence in SeqIO.parse(filename, handle):
        if sequence.id in selected_list:
            fasta_list.append(sequence)
    return fasta_list

def fastaList2Dic(list_fasta):
    """
    given a list of fasta, return a dic of fasta
    """
    dic_fasta = {}
    for sequence in list_fasta:
        dic_fasta[sequence.id] = sequence
    return dic_fasta

def open_fasta_to_dic(filename, handle = "fasta"):
    """
    given a file_name, open it with SeqIO, and make a list of fasta seqs
    """
    return SeqIO.to_dict(SeqIO.parse(filename, handle))

            

def fasta_sort_by_length(myfasta):
    """
    give a list of SeqIO input, return a list based on the length of seq
    """
    fastalen=fasta_length(myfasta)
    fastalength=fastalen.values()
    fastalength=list(set(fastalength))
    fastalength.sort()
    fastalendic={}
    fasta_sort=[]
    for item in fastalength:
        fastalendic[item]=[]
    for fasta in myfasta:
        fastalendic[len(fasta.seq)].append(fasta)
    #print(fastalendic)
    for item in fastalength:
        fasta_sort+=fastalendic[item]
    return fasta_sort
        
def fasta_unique_seq(myfasta):
    """
    give a list of SeqIO input, return a dictionary with ids of proteins with 
    the same sequence, and the amino acid sequence
    use the amino acid sequence as keys
    work for small list
    for big list, it might be slow
    use another function
    """
    uniquefasta={}
    numberseq=len(myfasta)
    for i in range(numberseq):
        aminoseq=str(myfasta[i].seq)
        if not aminoseq in uniquefasta:
            uniquefasta[aminoseq]=[myfasta[i].id]
        else:
            uniquefasta[aminoseq].append(myfasta[i].id)
    return uniquefasta



def fasta_uni_keepone(myfasta):
    """
    for sequence in myfasta, if more than one fasta with the same sequence, 
    keep only one of them. Return a list
    """
    fasta_seq_dic={}
    for sequence in myfasta:
        fasta_seq_dic[str(sequence.seq)]=sequence
    return list(fasta_seq_dic.values())

                    
def fasta_within(seq1,seq2list,nleft=0,nright=0):
    """
    seq1 and seq2list in SeqIO format
    nleft: number of amino acid /nucleotide neglect of seq1
    nright: number of amino acid /nucleotide neglect of seq1
    if seq1 without nleft and nright inside seq2list members, return True
    """
    if nright == 0:
        seq1short=seq1.seq[nleft:]
    else:
        seq1short=seq1.seq[nleft:-nright]
    
    for seq2 in seq2list:
        if seq1short in seq2.seq:
            return True
    return False

def fasta_within_with_description(seq1,seq2list,nleft=0,nright=0):
    """
    seq1 and seq2list in SeqIO format
    nleft: number of amino acid /nucleotide neglect of seq1
    nright: number of amino acid /nucleotide neglect of seq1
    if seq1 without nleft and nright inside seq2list members, return True
    return a string with description of seq1 in wich sequence of seq2list
    """
    if nright == 0:
        seq1short=seq1.seq[nleft:]
    else:
        seq1short=seq1.seq[nleft:-nright]
    
    for seq2 in seq2list:
        if seq1short in seq2.seq:
            return True, seq1.id+" in " + seq2.id
    return False, seq1.id+" not exist in the list "

def fasta_within_seq_big_with_description(myfasta,nleft=0,nright=0):
    """
    myfasta is a list of sequences of SeqIO.parse
    nleft: number of amino acid /nucleotide neglect of seq1
    nright: number of amino acid /nucleotide neglect of seq1
    return a list after removing sequence that is inside the others.
    default setting, nleft = nright=0
    if only provide nleft, than nright is considered the same as nleft
    also return a list of string with description
    """
    myfasta_sort = fasta_sort_by_length(myfasta)
    listlen=len(myfasta)
    myfasta_nowithin=[]
    describes=[]
    if nright == 0:
        nright = nleft
    for dummy_i in range(listlen-1):
        within, describe_within = fasta_within_with_description(myfasta_sort[dummy_i],myfasta_sort[dummy_i+1:],nleft,nright)
        if not within:
            myfasta_nowithin.append(myfasta_sort[dummy_i])
        else:
            describes.append(describe_within)
    myfasta_nowithin.append(myfasta_sort[-1])
    return myfasta_nowithin,describes


def fasta_within_seq_big(myfasta,nleft=0,nright=0):
    """
    myfasta is a list of sequences of SeqIO.parse
    nleft: number of amino acid /nucleotide neglect of seq1
    nright: number of amino acid /nucleotide neglect of seq1
    return a list after removing sequence that is inside the others.
    default setting, nleft = nright=0
    """
    myfasta_sort = fasta_sort_by_length(myfasta)
    listlen=len(myfasta)
    myfasta_nowithin=[]
    for dummy_i in range(listlen-1):
        if not fasta_within(myfasta_sort[dummy_i],myfasta_sort[dummy_i+1:],nleft,nright):
            myfasta_nowithin.append(myfasta_sort[dummy_i])
    myfasta_nowithin.append(myfasta_sort[-1])
    return myfasta_nowithin
    
def seq_extract(seqs):
    """
    read file "list.txt"
    extract sequence from seqs
    output to "seq_selected.txt", fasta format
    """
    #from Bio import SeqIO
    filetemp=open("list.txt","r")
    fo=open("seq_selected.txt","w")
    idstemp=list(filetemp)
    ids=[]
    for idtemp in idstemp:
        ids.append(idtemp[:-1])
    for seq in seqs:
        if seq.id in ids:
            SeqIO.write(seq,fo,"fasta")
    fo.close()
    filetemp.close()
    print("Done!")
    return


def saveFastaListToFile(mylist,myfile,myformat="fasta"):
    """mylist is a list of fasta sequences
    fo is the file name
    myformat is the saving format
    creating a file myfile, and save it to the working directory
    myformat is the same as those used in Biopython, SeqIO
    """
    fo=open(myfile,"w")
    for my_element in mylist:
        fo.write('>'+my_element.description+'\n'+str(my_element.seq)+'\n')
    fo.close()
    print("done!")    

def fasta_within_seq_big_fast(myfasta,nleft=0,nright=0):
    """
    myfasta is a list of sequences of SeqIO.parse
    nleft: number of amino acid /nucleotide neglect of seq1
    nright: number of amino acid /nucleotide neglect of seq1
    return a list after removing sequence that is inside the others.
    default setting, nleft = nright=0
    if only provide nleft, than nright is considered the same as nleft
    the same as fasta_within_seq_big, but this program can runs faster
    """
    myfasta_sort = fasta_sort_by_length(myfasta)
    listlen=len(myfasta)
    if listlen==1:
        return myfasta
    elif listlen==2:
        if myfasta_sort[0].seq in myfasta_sort[1].seq:
            return [myfasta_sort[1]]
        else:
            return myfasta_sort
    else:
         myfasta_nowithin=[]
    for dummy_i in range(listlen-1):
        if not fasta_within(myfasta_sort[dummy_i],myfasta_sort[dummy_i+1:],nleft,nright):
            myfasta_nowithin.append(myfasta_sort[dummy_i])
    myfasta_nowithin.append(myfasta_sort[-1])
    return myfasta_nowithin

def fasta_within_seq_big_faster(myfasta, nleft = 0, nright = 0):
    """
     myfasta is a list of sequences of SeqIO.parse
    nleft: number of amino acid /nucleotide neglect of seq1
    nright: number of amino acid /nucleotide neglect of seq1
    return a list after removing sequence that is inside the others.
    default setting, nleft = nright=0
    the same as fasta_within_seq_big, but this program can runs faster
    """
    myfasta_sort = fasta_sort_by_length(myfasta)
    listlen=len(myfasta)
    if listlen==1:
        return myfasta
    else:
        myfasta_nowithin = []
    allseqs = ""
    for fasta in myfasta_sort:
        allseqs += str(fasta.seq) + "\t"
    seqright = allseqs
    for dummy_i in range(listlen-1):
        (seqleft, seqright) = seqright.split("\t", 1)
        if nright >0:
            seqleftn = seqleft[nleft:-nright]
        if nright ==0:
            seqleftn = seqleft[nleft:]
        if seqright.count(seqleftn) == 0:
            myfasta_nowithin.append(myfasta_sort[dummy_i])
    myfasta_nowithin.append(myfasta_sort[-1])
    print('before',len(myfasta),' after',len(myfasta_nowithin))
    return myfasta_nowithin


def errorMatch(seq1, seq2, errors=2):
    """
    Given seq1 and seq2, len(seq1) <= len(seq2), return whether they match each other with allowed error.
    """
    if len(seq1) > len(seq2):
        return False
    step = len(seq1)//(errors+1)
    if step == 0:
        return True
    if errors == 0:
        return seq1 in seq2
    parts = [seq1[i:i+step] for i in range(0,len(seq1),step)] #separate seq1 to error+1 parts
    if len(parts[-1]) < step:
        parts[-2] = parts[-2]+parts[-1]
        parts.pop()
    similar = False
    sameslist =[]
    for i in range(errors+1):
        findsame = seq2.find(parts[i])
        if  findsame >= 0:
            similar = True
            sameslist.append((i,findsame))
    if similar == False:
        return False
    for i,j in sameslist:
        if j-step*(i)>=0 and j-step*(i)+len(seq1) <= len(seq2):
            seq2n = seq2[j-step*i:j-step*i +len(seq1)]
            missmatched = 0
            for k in range(len(seq1)):
                if seq1[k] != seq2n[k]:
                    missmatched += 1
                if missmatched > errors:
                    break
            if missmatched <= errors:
                return True
    return False

#errorMatch("ABCDEFG","ACCDEF")
#errorMatch("AB","ACCDGFF")
#errorMatch("BBCDEFG","ACCDEFG",1)
#errorMatch("ACCDEEGH","ABCDACCDEFGH",1)

        
    
def fasta_within_seq_big_withError(myfasta, error_rate = 0.02,kmerlen = 6):
    """
    myfasta is a list of SeqIO elements
    if a sequence is part of the other, with error_rate allowed, then remove this sequence.
    return a list of non-redundant SeqIO fasta
    """
    myfasta = fasta_uni_keepone(myfasta)
    import time
    time1 = time.time()
    dickmernum = {} #kmer dic, kmer with its seqs
    for dummyi in range(len(myfasta)):
        seq = str(myfasta[dummyi].seq)
        for i in range(len(seq)+1-kmerlen):
            kmernum = seq[i:i+kmerlen]
            if kmernum not in dickmernum:
                dickmernum[kmernum] = set()
            dickmernum[kmernum].add(dummyi)
    print(time.time()-time1)
    
    time1 = time.time() #change values of dickmernum to list
    for kmernum in dickmernum:
        dickmernum[kmernum] = list(dickmernum[kmernum])
    print(time.time()-time1)
    
    time1 = time.time()
    toremove = set()
    from collections import Counter
    for num1 in range(len(myfasta)):
        seq1 = str(myfasta[num1].seq)
        seq1kmers = set() # all kmernum, here is kmer5 in seq1
        for i in range(len(seq1)+1-kmerlen):
            seq1kmers.add(seq1[i:i+kmerlen])
    #    print(time.time()-time1)
        seq1targets = []
        for kmernum in seq1kmers:
            seq1targets += dickmernum[kmernum]
        seq1targets = Counter(seq1targets) # count the number of common kmers for each targets
        seq1targets = seq1targets.most_common() # sort the targets based on the number of commn kmers
    #    print(time.time()-time1)
        errors = int(len(seq1)*error_rate)
        for seq2id, seq2_counts in seq1targets:
            if seq2id != num1:
                seq2 = str(myfasta[seq2id].seq)
                if seq2_counts >= len(seq1kmers)-errors*kmerlen:
                    if len(seq1) <= len(seq2):
                        if seq2id not in toremove:
                            if seq2_counts >=10:
                                if errorMatch(seq1,seq2,errors):
                                    toremove.add(num1)
                                    break
    
    print(time.time()-time1)
    print(len(toremove))
    nonredunfasta =[]
    for i in range(len(myfasta)):
        if i not in toremove:
            nonredunfasta.append(myfasta[i])
    return nonredunfasta
    


