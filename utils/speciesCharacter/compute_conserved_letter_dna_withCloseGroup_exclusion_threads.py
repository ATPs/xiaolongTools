# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 16:04:08 2019

@author: ATPs
"""

#fa = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/samples_exons_for_paper.fasta'
#fn = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/group2_samplesID'
#fclose = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/all_samples_except_group2_and_outgroup'
#fbad = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/no_bad_samples'
#fout = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/outgroups'
#fexclude = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/zz'
#frange = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/fasta_dna_range'
#frename = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/pyrr_filtered_ordered_renamed.trelist'
#gapInGroup = 1.0 #percent of letters in ingroup
#gapInClose = 0.9 #percent of letters in closegroup
#gapOutGroup = 0.7 #percent of letters in the rest
#consOutGroup = 0.7 #the conservation of the most dominated letter in the outgroup




import os
from Bio import SeqIO
from collections import Counter
from multiprocessing import Pool
gapChars = set(['-', 'N'])
import numpy as np
import pickle
import uuid
import pandas as pd

file_shared_variables_in_memory = '/dev/shm/' + str(uuid.uuid4())
DO_NOT_KEEP_N_IN_TARGET = True


def read_nameDict(filename):
    '''
    filename stores the naming of sample_ID to annotated_ID
    each line looks like
    7993    7993_Celaenorrhinus syllius_Ecuador_1-3-Jan-2002
    return dict with k and v:
    7993    7993_Celaenorrhinus_syllius_Ecuador_1-3-Jan-2002
    '''
    dc = {}
    for line in open(filename,'r'):
        es = line.strip().split()
        dc[es[0]] = '_'.join(es[1:])
    return dc


def get_genus(adict):
    '''
    value is like 7993_Celaenorrhinus_syllius_Ecuador_1-3-Jan-2002
    return dict with Celaenorrhinus as value
    if not possible, use 'NA' as genus
    '''
    bdict = {}
    for k,v in adict.items():
        es = v.split('_')
        if len(es) > 1:
            genus = es[1]
        else:
            genus = 'NA'
        bdict[k] = genus
    return bdict


def read_IDs(filename):
    lines = open(filename).readlines()
    lines =[line.strip() for line in lines if line.strip() != '']#remove empty lines and \n symbol
    if len(lines) == 0:
        return set([])
    aset = set([line.strip().split()[0].split('_')[0] for line in lines])
    return aset


def collapse_seq_forGenus(seqlist):
    #this sequence is intentionally used for checking if the close species has a gap in the position
    #so when collapsing for genus, '-' is gap and 'N' is can not sure
    '''
    for each site, return '-' if with only gap chars '-' or 'N'
    return base if all are the same
    return N in other cases
    '''
    new = []
    for i in range(len(seqlist[0])):
        chars = [seq[i] for seq in seqlist if seq[i] not in gapChars]
        if len(chars) == 0:
            new.append('-')
        else:
            check = set(chars)
            if len(check) == 1:
                new.append(chars.pop())
            else:
                new.append('N')
    return ''.join(new)


def collapse_species(alist, gdict, badIDs):
    '''
    for each genus, return a "combined" sequence
    '''
    adict = {}
    for name, seq in alist:
        genus = gdict[name]
        try:
            adict[genus].append((name, seq))
        except KeyError:
            adict[genus] = [(name, seq)]

    rlist = []
    for name in adict:
        seqlist = adict[name]
        if len(seqlist) == 1:
            spID, seq = seqlist[0]
            if spID not in badIDs:
                rlist.append((name, seq))
        else:
            seqlist = [each[1] for each in seqlist]
            seq = collapse_seq_forGenus(seqlist)
            rlist.append((name, seq))
    return rlist


def read_range(fn):
    '''
    fn is a file of exon location of genes
    each line looks like
        #lac1000.1.e2	0	69
    for each ID, it has to include '.e' to indicate an exon
    '''
    df = pd.read_csv(fn,sep='\t',header=None)
    df.columns = ['exon_id','start','end']
    df['gene_id'] = df['exon_id'].apply(lambda x:x.split('.e')[0])
    df['exon_N'] = df['exon_id'].apply(lambda x:int(x.split('.e')[1]))
    df = df.sort_values(by = ['gene_id','exon_N'])
    dc = {} #gene_id as key, a list of start and end values of exons as value
    for _n, row in df.iterrows():
        if row['gene_id'] not in dc:
            dc[row['gene_id']] = []
        dc[row['gene_id']].append([row['start'],row['end']])
    
    return dc


def detect_groupping_characters(seqlist, otherlist, segClose, gene, badIDs, genus_Dict, gapInGroup = 1.0, gapInClose = 0.9, gapOutGroup = 0.7, consOutGroup = 0.7):
#    global genus_dict
    taken_chars = []
    #print seqlist
    seqlist_byGenus = collapse_species(seqlist, genus_Dict, badIDs)
    for index in range(len(seqlist[0][1])):
        #firstly take only the good samples and check gap and conservation
        chars = [seq[index] for name, seq in seqlist_byGenus if name not in badIDs]
        check_chars = [char for char in chars if char != '-'] #because gap has been collapsed, so N is meaningful here
        #if len(check_chars) == len(chars):
        Ningroup = len(chars) - len(check_chars)
        chars_target_sample = [seq[index] for name, seq in seqlist if name not in badIDs]
        Ningroup_by_sample = len([char for char in chars_target_sample if char == '-'])
        if len(check_chars) >= gapInGroup * len(chars):
            if len(set(check_chars)) == 1:
                goodchar = check_chars.pop()
                if DO_NOT_KEEP_N_IN_TARGET:
                    if goodchar == 'N':
                        continue

                #then it needs to check if the goodchar has confliction in samples
                isConsistent = True
                for name, seq in seqlist:
                    if name in badIDs:
                        char = seq[index]
                        if char not in gapChars:
                            if char != goodchar:
                                isConsistent = False
                if not isConsistent:
                    continue

                #the segClose is by genus
                #currently allow poor samples
                otherClose = set([seq[index] for name, seq in segClose])
                check_close = [char for char in otherClose if char != '-'] #close seq is also gap collpased, N is meaningful
                if len(check_close) < gapInClose * len(otherClose):
                    #didn't allow it to be gap
                    continue

                #the letter need to be different in other clade
                otherChars = [seq[index] for name, seq in otherlist if name not in badIDs]
                otherTotalChars = [seq[index] for name, seq in otherlist]

                check_others = [char for char in otherChars if char not in gapChars]
                #print 'consistent conserved letter found for %s %s %s' % (index, ','.join(check_chars), ','.join(otherChars))

                if len(check_others) > gapOutGroup * len(otherChars):
                    if goodchar in otherTotalChars:
                    #if check_others.count(goodchar) > 2:
                        continue
                    count_dict = Counter(check_others)
                    maxChar = max(count_dict, key=lambda x: count_dict[x])
                    if count_dict[maxChar] < consOutGroup * len(check_others):
                        continue
                    print('consistent conserved letter found for %s %s %s' % (index, ','.join(check_chars), ','.join(otherChars)))
                    fraction = check_others.count(goodchar) / float(len(check_others))
                    if fraction < 0.05:
                        #consP = count_dict[maxChar]/float(len(check_others))
                        consP = len(check_others) - count_dict[maxChar]
                        #gapP = len(check_others)/ float(len(otherChars))
                        gapP = len(otherChars) - len(check_others)
                        #taken_chars.append((index, goodchar, maxChar, Ningroup, consP, gapP))
                        taken_chars.append((index, goodchar, maxChar, Ningroup_by_sample, consP, gapP))

    #reformat
    if len(taken_chars) == 0:
        return 0, None
    else:
        N = 0
        info = []
        for index, char, other, Ningroup, consP, gapP in taken_chars:
            N += 1
            info.append(['%s%s%s' % (other, index+1, char), Ningroup, consP, gapP, gene, N])

        return N, info

def getSeq(seq, genelocs):
    '''
    seq is a sequence
    gene locs is a list of start and stop sites in seq
    return the sequence
    '''
    geneseq = []
    for start, end in genelocs:
        geneseq.append(seq[start:end])
    return ''.join(geneseq)

def processGene(ls_gene_info):
    ls_result = []
    with open(file_shared_variables_in_memory,'rb') as f:
        seqlist, otherlist, closeSeqList, badIDs, genus_Dict, gapInGroup, gapInClose, gapOutGroup, consOutGroup = pickle.load(f)
        
    for gene_info in ls_gene_info:
        gene,genelocs = gene_info
        print(gene, genelocs)
        print(gene,sum(e[1] for e in genelocs) - sum(e[0] for e in genelocs))
        seglist = [(name, getSeq(seq, genelocs)) for name, seq in seqlist]
        segOther = [(name,getSeq(seq, genelocs)) for name, seq in otherlist]
        segClose = [(name, getSeq(seq, genelocs)) for name, seq in closeSeqList]
        N, char_info = detect_groupping_characters(seglist, segOther, segClose, gene, badIDs, genus_Dict=genus_Dict, gapInGroup = gapInGroup, gapInClose = gapInClose, gapOutGroup = gapOutGroup, consOutGroup = consOutGroup)
        ls_result.append([N,char_info])
    return ls_result


def computeConservedLetterDNAwithCloseGroupExclusion(fa, fn, fclose, fbad, fout, fexclude, frange, frename, gapInGroup=1.0, gapInClose=0.9, gapOutGroup=0.7, consOutGroup=0.7, outfilename = None, threads = 16):
    '''
    fa: fasta alignment of all samples
    fn: a file with the sample_ids of target samples to check
    fclose: a file with sample_ids of close related species.
    fbad: a file with sample_ids of bad samples to be excluded from the analysis of the ingroups
    fout: a file with sample_ids of outgroup samples
    fexclude: a file with samples_ids of samples to be exculded. based on wenlin's code, fexclude is merged is fout
    frange: a file with the locatio of exons in the fasta alignment.
            each line with three elements: exon_id, start_position, end_position.
            exon_id must be in this format: XXXXXXX.eN, '.e' is essential
    frename: a file of sampl_id and long description for that sample
    gapInGroup: percent of letters in ingroup. a float number to filter the gap ratio of InGroup samples. default 1.0
    gapInClose: percent of letters in closegroup. a float number to filter the gap ratio of InClose samples. default 0.9
    gapOutGroup: percent of letters in the rest. a float number to filter the gap ratio of OutGroup samples. default 0.7
    consOutGroup: the conservation of the most dominated letter in the outgroup. a float number. default 0.7 
    '''
    #nameDict, sample_id to long names with description
    print('read in nameDict')
    nameDict = read_nameDict(filename=frename)
    print('nameDict:')
    for k,v in nameDict.items():
        print(k,v)
    
    #genus_Dict, sample_id to species genus. In default, the input nameDict have value like 7993_Celaenorrhinus_syllius_Ecuador_1-3-Jan-2002, if split with '_', the second element is genus
    genus_Dict = get_genus(nameDict)
    
    #read in sampleIDs in different groups
    takenIDs = read_IDs(fn)
    if fbad is None:
        badIDs = set([])
    else:
        badIDs = read_IDs(fbad)
    closeIDs = read_IDs(fclose)
    outIDs_1 = read_IDs(fout)
    if fexclude is None:
        excludeIDs = set([])
    else:
        excludeIDs = read_IDs(fexclude)
    outIDs = outIDs_1 | excludeIDs
    print('\n\n    read in input ID files')
    print('for these file, each line includes a sample_ID')
    print('number of takenIDs:', len(takenIDs))
    print('number of badIDs', len(badIDs))
    print('number of closeIDs', len(closeIDs))
    print('number of input outIDs',len(outIDs_1))
    print('number of excludeIDs', len(excludeIDs))
    print('excludeIDs and input outIDs merged, total outIDs', len(outIDs))
    
    #read in alignment fasta file
    ls_seqs = list(SeqIO.parse(fa,'fasta'))
    oSeqDict = {e.description:str(e.seq) for e in ls_seqs}
    seqDict = {name: oSeqDict[name] for name in oSeqDict if name.split('_')[0] not in outIDs}
    print('filtering for outgroups %s %s' % (len(oSeqDict), len(seqDict)))
    
    
    seqlist = [(name.split('_')[0], seqDict[name]) for name in seqDict if name.split('_')[0] in takenIDs]
    otherlist = [(name.split('_')[0], seqDict[name]) for name in seqDict if name.split('_')[0] not in takenIDs]
    closeSeqList = [(name.split('_')[0], seqDict[name]) for name in seqDict if name.split('_')[0] in closeIDs]
    #TODO: not counting bad samples in close group?
    closeSeqList = collapse_species(closeSeqList, genus_Dict, badIDs) #a list with (genus, combined_sequence). samples in badIDs will be excluded from the analysis
    
    
    print('number of sequence:', len(seqlist))
    print('inSeqNames:', ','.join([each[0] for each in seqlist]))
    print('number of other:', len(otherlist))
    if len(seqlist) != len(takenIDs):
        print('''Error! ID number in takenIDs does not agree with sequence number! Likely that some IDs are missing in the input fasta sequences''')
        print('missingIDs:')
        print('\n'.join(takenIDs - set([name.split('_')[0] for name in seqDict])))
    
    rangeDict = read_range(frange)
    dc_genelen = {}
    for gene, genelocs in rangeDict.items():
        dc_genelen[gene] = sum(e[1] for e in genelocs) - sum(e[0] for e in genelocs)
    
    if outfilename is None:
        dn = '%s_%s_correctIndex.cc' % (os.path.basename(fa), os.path.basename(fn))#where to store the output file
    else:
        dn = outfilename
    dp = open(dn, 'w')
    header = 'gene\tmutation_in_gene\tchar\tgapInGroup\tgapInRest\tnonConserveLetterInRest\n'
    dp.write(header)
    ls_gene_frag = list(rangeDict.items())
    print('numbers of genes to be processed:', len(ls_gene_frag))
    
    n = len(ls_gene_frag)
    step = int(np.ceil(n/threads))
    lsls_gene_frag = [ls_gene_frag[i:i+step] for i in range(0,n,step)]
    with open(file_shared_variables_in_memory,'wb') as f:
        pickle.dump([seqlist, otherlist, closeSeqList, badIDs, genus_Dict, gapInGroup, gapInClose, gapOutGroup, consOutGroup],f)
    del seqlist, otherlist, closeSeqList, oSeqDict, seqDict, ls_seqs
    
    
    pool = Pool(threads)
    lsls_results = pool.map(processGene, lsls_gene_frag)
    pool.close()
    ls_results = []
    os.remove(file_shared_variables_in_memory)
    for _ls_results in lsls_results:
        ls_results = ls_results + _ls_results
#        N, char_info = detect_groupping_characters(seglist, segOther, segClose, gene, badIDs, genus_Dict=genus_Dict, gapInGroup = gapInGroup, gapInClose = gapInClose, gapOutGroup = gapOutGroup, consOutGroup = consOutGroup)
    
    allInfo = []
    for N, char_info in ls_results:
        if N != 0:
            #print '\n'.join([each[1] for each in seglist])
            #print '\n'.join([each[1] for each in segOther])
            #print('++++++++++++++++++++')
            for each in char_info:
                each[-1] = N
            allInfo += char_info
    
    char_info = sorted(allInfo, key=lambda x: (x[1], -x[-1]))
    for char, Ningroup, consP, gapP, gene, N in char_info:
        line = '%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, N, char, Ningroup, gapP, consP)
        dp.write(line)
    dp.close()
    

description = '''
a function to find the conserved letter in sequence alignment. The inputs are.
    fa: fasta alignment of all samples
    fn: a file with the sample_ids of target samples to check
    fclose: a file with sample_ids of close related species. Typically it is those samples other than the target outgroup samples
    fbad: a file with sample_ids of bad samples to be excluded from the analysis of the ingroups
    fout: a file with sample_ids of outgroup samples
    fexclude: a file with samples_ids of samples to be exculded. based on wenlin's code, fexclude is merged with fout
    frange: a file with the locatio of exons in the fasta alignment.
            each line with three elements: exon_id, start_position, end_position.
            exon_id must be in this format: XXXXXXX.eN, '.e' is essential
    frename: a file of sampl_id and long description for that sample
    gapInGroup: percent of letters in ingroup. a float number to filter the gap ratio of InGroup samples. default 1.0
    gapInClose: percent of letters in closegroup. a float number to filter the gap ratio of InClose samples. default 0.9
    gapOutGroup: percent of letters in the rest. a float number to filter the gap ratio of OutGroup samples. default 0.7
    consOutGroup: the conservation of the most dominated letter in the outgroup. a float number. default 0.7 

below is how it works:
    first read in nameDict
    get sample_id to species genus based on nameDict. In default, the input nameDict have value like 7993_Celaenorrhinus_syllius_Ecuador_1-3-Jan-2002, if split with '_', the second element is genus
    read in sampleIDs in different input groups: fn, fclose, fbad, fout, fexclude
        currently, fout is combined with fexclude
    read in the sequence alignment
    split the alignment to 
        seqlist: sequences in fn, target group
        otherlist: sequences not in fn
        closeSeqList: sequences in fclose
            closeSeqList is further processed. 
            Sequences in fbad will be excluded from the analysis.
            for sequences from different genus, a concensus sequence will be reported. 
                the concensus sequences is produced as: for each site, if all sites of sequences in that genus are gap ('N' or '-'), use '-'. If there is only one kind of base of 'ATCG', use that base, otherwise, use 'N' to indicate not sure about the bases
    read in frange, the location of exons in fa, the alignment file. 
        return a dictionary with gene_id as key, [start, end] as value
        Note, for frage, each line with three elements: exon_id, start_position, end_position.
        exon_id must be in this format: XXXXXXX.eN, '.e' is essential. program with split exon_id by '.e' and get the location of whole gene. XXXXXXX is considered as gene_id
    process the data:
        for each gene in frange, get seglist (seqlist), segOther (otherlist), segClose(closeSeqList)
        then detect_groupping_characters
            for each gene, seqlist is also collapsed first by genus, excluding those badIDs
            for each site in the genes, first remove sites that is "-". Note, site with "N" is kept. It means "uncertain" instead of a gap here.
            it is considered a good site (goodChar) if:
                ratio of non-gap chars in "ingroup" (target group) >= gapInGroup
                there is only one kind of base in "ingroup"
                in "otherClose" group, the ratio of non-gap chars >= gapInClose
                in "otherlist", any samples not in the target group, ratio of non-gap chars >= gapOutGroup
                the base is not in outgroup bases
                the dominant base in outgroup is >= consOutGroup * number of non-gap bases in outgroup
                the fraction of target char is less than 0.05 of chars in Note: seems a redundant logic.
                return location of goodChar, goodChar base, most abundant char is otherlist, number of gaps in ingroup or other group, and other information


'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    print(description)
    parser.add_argument('-a','--fa', help = 'fasta alignment of all samples', required=True)
    parser.add_argument('-n','--fn', help = 'a file with the sample_ids of target samples to check', required=True)
    parser.add_argument('-c','--fclose', help = 'a file with sample_ids of close related species', required=True)
    parser.add_argument('-b','--fbad', help = 'a file with sample_ids of bad samples to be excluded from the analysis of the ingroups. default None, no samples are bad', default=None)
    parser.add_argument('-o','--fout', help = 'a file with sample_ids of outgroup samples', required=True)
    parser.add_argument('-e','--fexclude', help = '''a file with samples_ids of samples to be exculded. based on wenlin's code, fexclude is merged with fout. default None, no samples to exclude''', default=None)
    parser.add_argument('-r','--frange', help = 'a file with the locatio of exons in the fasta alignment', required=True)
    parser.add_argument('-R','--frename', help = 'a file of sampl_id and long description for that sample.', required=True)
    parser.add_argument('-I','--gapInGroup', help = 'percent of letters in ingroup. a float number to filter the gap ratio of InGroup samples. default 1.0', type=float, default=1.0)
    parser.add_argument('-C','--gapInClose', help = 'percent of letters in closegroup. a float number to filter the gap ratio of InClose samples. default 0.9', type=float, default=0.9)
    parser.add_argument('-O','--gapOutGroup', help = 'percent of letters in the rest. a float number to filter the gap ratio of OutGroup samples. default 0.7', type=float, default=0.7)
    parser.add_argument('-S','--consOutGroup', help = 'the conservation of the most dominated letter in the outgroup. a float number. default 0.7', type=float, default=0.7)
    parser.add_argument('-z','--outfilename', help = '''where to store the result outfile. default None, which will be (fa)_(fn)_correctIndex.cc''', default=None)
    parser.add_argument('-T','--threads', help = '''number of CPUs to use. default 16''', default=16, type=int)
    f = parser.parse_args()
    computeConservedLetterDNAwithCloseGroupExclusion(fa=f.fa, fn=f.fn, fclose=f.fclose, fbad=f.fbad, fout=f.fout, fexclude=f.fexclude, frange=f.frange, frename=f.frename, gapInGroup=f.gapInGroup, gapInClose=f.gapInClose, gapOutGroup=f.gapOutGroup, consOutGroup=f.consOutGroup, outfilename=f.outfilename,threads=f.threads)



