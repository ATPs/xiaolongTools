from Bio.SeqIO import parse

import io
import os
import sys

#used by read_fasta
import cmn

#used by get_id
import re

#a class to store data.
#aln is a list of SeqRecord Object from SeqIO that has attributes 'id,name,description,seq'
#data_type is the format of the data
#usefull function: fasta,tab;easy_index;cut_by_range;

#used by alignment.conseunsus()
from collections import Counter
#used in diff to determine the type
from Bio.SeqRecord import SeqRecord
#used by random_sample
from random import randint


class alignment(object):
    def __init__(self, input_aln, data_type='fasta'):
        if os.path.isfile(input_aln):
            info = self._read_file_(input_aln)
        else:
            info = input_aln
        self.data_type = data_type
        try:
            records = parse(io.StringIO(info), data_type)
            aln = [entry for entry in records]
        except:
            aln = []
        if aln == []:
            print('\nError: cannot read the input. May be format problem.\nThe format labels are defined in http://biopython.org/wiki/AlignIO\n', file=sys.stderr)
            #sys.exit()
        self.aln_list = aln

    def _read_file_(self, fn):
        f = open(fn)
        a = f.read()
        f.close()
        return a

    def isfasta(info):
        lines = info.split('\n')
        found = False
        line = 'go'
        k = len(lines)
        i = 0
        while (not found and i < k):
            line = lines[i]
            if line.startswith('>'):
                found = True
            i += 1
        return found

    def cut_by_range(self, cut_range, remove_all_gap=True):
        #cut the alignment by the input range
        #range should be like [[4,6],[4,5]] which is used in the list object

        def cutting(seq, cut_range):
            #example: [10,20] will cut 11 residues from 10 to 20 (included)
            new = seq[0:0]
            for cut in cut_range:
                if len(cut) == 1:
                    new += seq[:cut[0]]
                else:
                    new += seq[:cut[0]] + seq[cut[1] + 1:]
            return new

        cut_aln_list = list(self.aln_list)
        for i, aln in enumerate(self.aln_list):
            cut_aln_list[i].seq = cutting(aln.seq, cut_range)

        #reomve all gap one
        new_list = []
        if remove_all_gap:
            for aln in cut_aln_list:
                if str(aln.seq).replace('-', '').strip() != '':
                    new_list += [aln]
            cut_aln_list = new_list
        #regenerate the alignment object
        alist = ['>%s\n%s\n' % (aln.description, aln.seq)
                 for aln in cut_aln_list]
        return alignment(''.join(alist))

    def fasta(self, number=None, aln=None):
        if aln!=None:
            return '>%s\n%s\n' % (aln.description.seq)
        elif number is None:
            return ['>%s\n%s' % (aln.description, aln.seq) for aln in self.aln_list]
        else:
            return ['>%s\n%s' % (aln.description, aln.seq)
                    for aln in self.aln_list[number]]

    def tab(self):
        return ['%s\t%s' % (aln.description, aln.seq) for aln in self.aln_list]

    def easy_index(self, col=None):
        ########################
        def make_head(length):
            line1 = ''
            line2 = ''
            for i in range(length):
                if i % 10 == 0:
                    a = str(i / 10)
                    line1 += a
                    number_space = 10 - len(a)
                    line1 += ' ' * number_space

                line2 += str(i % 10)
            return [line1, line2]
        ###########################

        seq_length = len(self.aln_list[0].seq)
        heads = make_head(seq_length)
        labels = [aln.name[:10] for aln in self.aln_list]
        seqs = [aln.seq for aln in self.aln_list]

        if col is None:
            flist = [' ' * 14 + i + '\n' for i in heads]
            flist += ['%s    %s\n' % (labels[i], seqs[i])
                      for i in range(len(seqs))]
        else:
            flist = []
            for step in range((seq_length / col) + 1):
                #print step
                flist += [' ' * 14 + line[step * col:(step + 1) * col] + '\n'
                          for line in heads]
                flist += ['%s    %s\n' % (labels[i], seqs[i][step * col:(step + 1) * col])
                          for i in range(len(seqs))]
                flist += ['\n']

        return ''.join(flist)

    def size(self):
        return [len(self.aln_list), len(self.aln_list[0].seq)]

    def regenerate_alignment(self, aln_list):
        alist = ['>%s\n%s\n' % (aln.description, aln.seq)
                 for aln in aln_list]
        return alignment(''.join(alist))

    def remove_by_gap(self, gap_fraction=0.5):
        record = []
        for aln in self.aln_list:
            seq = str(aln.seq)
            gaps = float(seq.count('-')) / len(seq)
            if gaps < gap_fraction:
                record += [aln]
        if record == []:
            print('no sequence left after cutting the gap...')
            return []

        return self.regenerate_alignment(record)

    def pairwise_identity(self,aln1=None,aln2=None,skip_gap=True):
        """
        returned value (x) range from 0 to 100. (x%)
        identity calculated by identical positions divided by the longest sequence
        skipping gap doesn't affect the total length
        """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #this is used for comparing two seq
        if aln1!=None and aln2!=None:
            #aln1,2 is the SeqIO object
            #if not return_position:
            if True:
                seq1=str(aln1.seq)
                seq2=str(aln2.seq)
                identity=0
                kk=len(seq1)
                if len(seq2)>kk:
                    kk=len(seq2)
                for ii in range(kk):
                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1
                    #else:
                    #    print seq1[ii],seq2[ii]
                #print seq1,'\n',seq2
                return round(identity*100.0/kk,1)

            #else:#return position
             #   seq1=str(aln1.seq)
              #  seq2=str(aln2.seq)
               # identity=0
                #kk=len(seq1)
                #diff_position=[]
                #for ii in xrange(kk):
                #    if seq1[ii]==seq2[ii]:
                 #       if not skip_gap:
                  #          identity+=1
                   #     elif seq1[ii]!='-' and seq2[ii]!='-':
                    #        identity+=1
                    #else:#seq is different (not include gap)
                     #   diff_position.append(ii)
            #return [kk-identity,diff_position]
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        adict={}
        k=len(self.aln_list)
        for i in range(k):
            for j in range(i+1,k):
                seq1=str(self.aln_list[i].seq)
                seq2=str(self.aln_list[j].seq)
                kk=len(seq1)
                identity=0

                if not skip_gap:
                    length=kk
                else:
                    length=0#if skip gap, only count the region both having char

                for ii in range(kk):
                    if skip_gap:
                        if seq1[ii]!='-' and seq2[ii]!='-':
                            length+=1

                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1

                adict['%s %s' % (i,j)]=str(round(identity*100.0/length,1))
        return adict

    #used once
    #calculate the identity to be integer, only record the frequency
    def pairwise_identity_stat(self,skip_gap=True):
        """
        calculate the pairwise identity for all sequence in a alignment
        the output file is 'idt number'. for example,
        '80 765' means there are 765 sequence pairs whose identity is 80%
        """
        adict={}
        k=len(self.aln_list)
        for i in range(k):
            for j in range(i+1,k):
                seq1=str(self.aln_list[i].seq)
                seq2=str(self.aln_list[j].seq)
                kk=len(seq1)
                identity=0

                if not skip_gap:
                    length=kk
                else:
                    length=0#if skip gap, only count the region both having char

                for ii in range(kk):
                    if skip_gap:
                        if seq1[ii]!='-' and seq2[ii]!='-':
                            length+=1

                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1
                key=int(identity*100/length)
                if key not in adict:
                    adict[key]=1
                else:
                    adict[key]+=1
                #adict['%s %s' % (i,j)]=str(round(identity*100.0/length,1))
        return adict

    def pairwise_difference(self,skip_gap=True,aln1=None, aln2=None,return_position=False):
        """
        used to calculate the NUMBER of DIFFerent chars in pairwise sequences
        if aln1 and aln2 indicated, only calcalute these two sequence,
        else, calculate all sequences in the my_aln object
        return_position is only validated in the 'aln1 and aln2' mode
        """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #this is used for comparing two seq
        if aln1!=None and aln2!=None:
            #aln1,2 is the SeqIO object
            if not return_position:
                seq1=str(aln1.seq)
                seq2=str(aln2.seq)
                identity=0
                kk=len(seq1)
                for ii in range(kk):
                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1
                    #else:
                    #    print seq1[ii],seq2[ii]
                #print seq1,'\n',seq2
                return kk-identity

            else:#return position
                seq1=str(aln1.seq)
                seq2=str(aln2.seq)
                identity=0
                kk=len(seq1)
                diff_position=[]
                for ii in range(kk):
                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1
                    else:#seq is different (not include gap)
                        diff_position.append(ii)
            return [kk-identity,diff_position]
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        adict={}
        k=len(self.aln_list)
        for i in range(k):
            for j in range(i+1,k):
                seq1=str(self.aln_list[i].seq)
                seq2=str(self.aln_list[j].seq)
                kk=len(seq1)
                identity=0

                if not skip_gap:
                    length=kk
                else:
                    length=0#if skip gap, only count the region both having char

                for ii in range(kk):
                    if skip_gap:
                        if seq1[ii]!='-' and seq2[ii]!='-':
                            length+=1

                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1
                    #else:
                    #    print seq1[ii],seq2[ii]
                #print seq1,'\n',seq2
                #print i,j,identity, length
                adict['%s %s' % (i,j)]=str(length-identity)
        return adict

    def pairwise_difference_stat(self,skip_gap=True):
        """
        similar to pairwise_identity_stat, but
        this function calculate the NUMBER of identical chars,
        but pairwise_identity_stat uses the PERCENTAGE of identity as key
        """
        adict={}
        k=len(self.aln_list)
        for i in range(k):
            for j in range(i+1,k):
                seq1=str(self.aln_list[i].seq)
                seq2=str(self.aln_list[j].seq)
                kk=len(seq1)
                identity=0

                if not skip_gap:
                    length=kk
                else:
                    length=0#if skip gap, only count the region both having char

                for ii in range(kk):
                    if skip_gap:
                        if seq1[ii]!='-' and seq2[ii]!='-':
                            length+=1

                    if seq1[ii]==seq2[ii]:
                        if not skip_gap:
                            identity+=1
                        elif seq1[ii]!='-' and seq2[ii]!='-':
                            identity+=1

                key=identity
                if key not in adict:
                    adict[key]=1
                else:
                    adict[key]+=1
        return adict


    def consensus(self,internal=False):
        #get the consensus residue for each position
        #index starting with 0
        #output as a dictionary
        adict={}
        M,N=self.size()
        alist=[]
        for i in range(N):
            tmp=[self.aln_list[x].seq[i] for x in range(M)]
            freq=Counter(tmp).most_common(1)[0]
            alist+=[freq[0]]
            adict[i]={'residue':freq[0],'frequency':float(freq[1])/M}
        if not internal:
            return adict
        else:#used by consensus_sequences, output the consensus sequence
            return adict,''.join(alist)

    #return the number of different residues
    def diff(self,a,b):
        if isinstance(a,SeqRecord):
            seq1=str(a.seq)
        elif isinstance(a,str):
            seq1=a

        if isinstance(b,SeqRecord):
            seq2=str(b.seq)
        elif isinstance(b,str):
            seq2=b

        seq=list(zip(seq1,seq2))
        count=0
        for pair in seq:
            if pair[0]!=pair[1]:
                count+=1
        return count

    def consensus_sequences(self,check_close=True):
        #check_close: get the sequence closest to the consensus sequence
        #if False, closest=''
        adict,seq=self.consensus(internal=True)
        if check_close:
            min_diff=9999
            for aln in self.aln_list:
                diff_res=self.diff(aln,seq)
                if diff_res<min_diff:
                    min_diff=diff_res
                    closest=str(aln.seq)
            return [seq,closest]
        else:
            return [seq,'']


    def unify(self):
        #remove duplication. for same sequence, only keep one.
        seq_set=set([])
        new_list=[]
        for aln in self.aln_list:
            seq=str(aln.seq)
            if seq not in seq_set:
                seq_set.add(seq)
                new_list.append(aln)
        self.aln_list=new_list
        return None

    def Meff(self,idt=80,aln=None,fast=True):
        """
        culculate the effective sequence number of the alignment;
        default cutoff is 80 (>80% identity)
        input aln is the same format as self.aln_list
        if False==True, use the fast_differ to do, such treat gap as the 21th type
        """
        if aln==None:
            thealn=self.aln_list
        else:#aln!=None
            thealn=aln

        M=0
        for aln1 in thealn:
            count=1
            #print 'new identity!\n'
            for aln2 in thealn:
                if aln2!=aln1:#already counted in the count
                    if fast==False:
                        identity=self.pairwise_identity(aln1=aln1,aln2=aln2)
                        #print identity
                        if identity>idt:
                            count+=1
                    else:
                        diff=self.fast_differ(aln1.seq,aln2.seq)
                        if (1-diff)*100>idt:
                            count+=1
            W=1.0/count
            M+=W
        return M

    def random_sample(self,aln=None,size=1):
        """
        randomly sample the sequences from the alignment
        the resulted set is of the size (default=1)
        return the list consisting of SeqRecord objects
        """
        if aln==None:
            thealn=self.aln_list
        else:
            thealn=aln

        M=len(thealn)
        if size>M:
            raise IndexError('input size larger than the size of alignment')
        final_set=set([])
        M=M-1
        while (len(final_set)<size):
            i=randint(0,M)
            final_set.add(thealn[i])
        return final_set

    def fast_differ(self,aln1,aln2):
        """
        aln1 and aln2 are two string-like object
        return the degree of difference (0--1)
        """
        #used by fast_pairwise_difference
        scipy_path='/usr6/local/lib/python2.7/site-packages/'
        if scipy_path not in sys.path:
            sys.path.append(scipy_path)
        from scipy import spatial
        not_identity=spatial.distance.hamming(aln1,aln2)
        return not_identity



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_fasta(fn):
    """
    read the fasta file, put it into SeqRecord Object
    """
    if True:
        if os.path.isfile(fn):
            info = cmn.txt_read(fn)
        else:
            info = fn
        if True:
            records = parse(io.StringIO(info), 'fasta')
            aln = [entry for entry in records]
        return aln

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def output_fasta(sth):
    """
    reverse version of read_fasta, put SeqRecord into string
    """
    if isinstance(sth,SeqRecord):
        return '>%s\n%s' %  (sth.description,sth.seq)
    elif isinstance(sth,list):
        alist=[]
        for i in sth:
            if isinstance(i,SeqRecord):
                alist.append('>%s\n%s' % (i.description,i.seq))
            else:
                print(sys.stderr, 'Only accept SeqRecord Object...Omitting record %s' % i)
        return '\n'.join(alist)

    else:
        print(sys.stderr, 'type %s is not acceptable...' % type(sth))
        return ''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_id(label,defline):
    m=re.compile('%s\|(.+?)\|' % label)
    return m.search(defline).group(1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
