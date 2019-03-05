# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 11:12:28 2019

@author: ATPs
"""

#fn ='/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/samples_exons_for_paper.fasta_group2_samplesID_correctIndex.cc'
#faln = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/samples_exons_for_paper.fasta'
#ftree= '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/Pyrrhopyginae_nuclear_for_paper.rescaled.tre'
#fbad = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/no_bad_samples'
#fnameTable = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/pyrr_filtered_ordered_renamed.trelist'
#frange = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/fasta_dna_range'
#f_protein_annotation = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/aly_table'



from collections import Counter
import ete3
import os
from Bio import SeqIO
import pandas as pd


def order_name_by_tree(names, tree):
    #midpoint = tree.get_midpoint_outgroup()
    #tree.set_outgroup(midpoint)

    names = [node.name for node in tree]
    new = []
    for name in names:
        name = name.strip("'").split('_')[0]
        new.append(name)
        #if name == '5727_dna':
        #    new.append('14063C11_dna')
    return new


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


def rename_names(names, badIDs, fnameTable):
    #fnameTable = 'allIDlist.renamed'
    nameDict = {}
    for line in open(fnameTable):
        if line.strip() == '':
            continue
        items = line.strip().split()
        nameDict[items[0]] = '_'.join(items[1:])

    new = []
    for name in names:
        name = name.split('_')[0]
        newName = nameDict[name]
        #if name in badIDs:
        if name in badIDs:
            newName = '[gapRM]%s' % newName
        new.append(newName)
    return new


def format_names(names, badIDs,fnameTable):
    #TODO: mark deleted names
    renamed_names = rename_names(names, badIDs, fnameTable)
    #print 'renamed_names', renamed_names
    maxLength = max(map(len, renamed_names))
    formater = '{:<%s}' % maxLength
    new = {}
    for name, newName in zip(names, renamed_names):
        new[name] = formater.format(newName)

    return new


def read_anno_table(fanno):
    nameDict = {}
    for line in open(fanno):
        items = line.strip().split('\t')
        nameDict[items[0]] = items[1]

    return nameDict


def get_alignment_range(chars):
    index = int(chars[1:-1]) - 1
    i = index - 10
    j = index + 11
    return i, j

def read_IDs(filename):
    lines = open(filename).readlines()
    lines =[line.strip() for line in lines if line.strip() != '']#remove empty lines and \n symbol
    if len(lines) == 0:
        return set([])
    aset = set([line.strip().split()[0].split('_')[0] for line in lines])
    return aset


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


def prepareRegionsAndAnno(fn, faln, ftree, fbad, fnameTable, frange, f_protein_annotation, f_ingroupNames = None):
    '''
    fn is the result of compute_conserved_letter_dna_withCloseGroup_exclusion.
    faln is the file of sequence alignment
    ftree is a tree file which will be used to re-order samples
    fbad is a file with bad sample_ids that should be excluded from the analysis
    fnameTable is a file of sampl_id and long description for that sample
    frange: a file with the locatio of exons in the fasta alignment.
        each line with three elements: exon_id, start_position, end_position.
        exon_id must be in this format: XXXXXXX.eN, '.e' is essential
    f_protein_annotation: a file with protein_id and annotation. typically a result of protein function annotation. the first column is protein_id and second column is function
    f_ingroupNames is a file with the target samples. if provided, output alignments with target samples first
    '''
    #read in fasta alignments
    fa = faln
    ls_seqs = list(SeqIO.parse(fa,'fasta'))
    seqDict = {e.description.split('_')[0]:str(e.seq) for e in ls_seqs} #only use sample_id as key
    
    #order names in seqDict by the tree file
    ordered_names = order_name_by_tree(seqDict, ete3.Tree(ftree))
    print('ordered_names', ordered_names)
    
    #read frange
    rangeDict = read_range(frange)
    dc_genelen = {}
    for gene, genelocs in rangeDict.items():
        dc_genelen[gene] = sum(e[1] for e in genelocs) - sum(e[0] for e in genelocs)
    
    badIDs=set([])
    printName_dict = format_names(ordered_names, badIDs, fnameTable)
    #print 'printName_dict', printName_dict
    for name in printName_dict:
        print('checkPrint', name, '"%s"' % printName_dict[name])
    
    geneAnno = read_anno_table(f_protein_annotation)
    
    dn = '%s.alignment' %fn
    dlist = '%s.anno' %fn
    lines = list(open(fn))
    lines = [e.strip() for e in lines]
    with open(dn, 'w') as dp, open(dlist, 'w') as dplist:
        dplist.write('%s\tannotation\n' % lines[0])
        for line in lines[1:]:
            print(line)
            items = line.strip().split()
            gene, N, chars, Nin, gapP, consP = items[:6]
            genelocs = rangeDict[gene]
            print(gene, dc_genelen[gene], genelocs)
            if True:
                anno = geneAnno[gene]
                dplist.write('%s\t"%s"\n' % (line, anno))
                dp.write('\n\n###alignment for %s %s (gapP: %s; SNPs: %s)###\n' % (gene, chars, gapP, N))
                ii, jj = get_alignment_range(chars)
                for name in ordered_names:
                    if name not in seqDict:
                        print('warning!!! error!!!', name, 'not found in the input sequence alignment')
                        continue
                    geneSeq = getSeq(seqDict[name], genelocs)
                    line_2 = '%s    %s\n' % (printName_dict[name], geneSeq[ii:jj])
                    dp.write(line_2)
    
    if f_ingroupNames is not None:
        print('target sample name provided, also output alignment with target group first')
        takenIDs = read_IDs(f_ingroupNames)
        takenIDs = list(takenIDs)
        new_ordered_names1 = [name for name in takenIDs if name in ordered_names]
        new_ordered_names = new_ordered_names1 + [name for name in ordered_names if name not in new_ordered_names1]
        dn = '%s.alignment.sorted' %fn
        lines = list(open(fn))
        lines = [e.strip() for e in lines]
        with open(dn, 'w') as dp:
            for line in lines[1:]:
                print(line)
                items = line.strip().split()
                gene, N, chars, Nin, gapP, consP = items[:6]
                genelocs = rangeDict[gene]
                print(gene, dc_genelen[gene], genelocs)
                if True:
                    anno = geneAnno[gene]
                    dp.write('\n\n###alignment for %s %s (gapP: %s; SNPs: %s)###\n' % (gene, chars, gapP, N))
                    ii, jj = get_alignment_range(chars)
                    for name in new_ordered_names:
                        if name not in seqDict:
                            print('warning!!! error!!!', name, 'not found in the input sequence alignment')
                            continue
                        geneSeq = getSeq(seqDict[name], genelocs)
                        line_2 = '%s    %s\n' % (printName_dict[name], geneSeq[ii:jj])
                        #print(line_2)
                        dp.write(line_2)
    
    

description = '''
fn is the result of compute_conserved_letter_dna_withCloseGroup_exclusion.
faln is the file of sequence alignment
ftree is a tree file which will be used to re-order samples
fbad is a file with bad sample_ids that should be excluded from the analysis
fnameTable is a file of sampl_id and long description for that sample
frange: a file with the location of exons in the fasta alignment.
    each line with three elements: exon_id, start_position, end_position.
    exon_id must be in this format: XXXXXXX.eN, '.e' is essential
f_protein_annotation: a file with protein_id and annotation. typically a result of protein function annotation. the first column is protein_id and second column is function, delimiter is tab sysmbol 
'''
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    print(description)
    parser.add_argument('-n','--fn', help = 'fn is the result of compute_conserved_letter_dna_withCloseGroup_exclusion', required=True)
    parser.add_argument('-a','--faln', help = 'the file of sequence alignment', required=True)
    parser.add_argument('-t','--ftree', help = 'a tree file which will be used to re-order samples', required=True)
    parser.add_argument('-b','--fbad', help = 'a file with bad sample_ids that should be excluded from the analysis, default None', default=None)
    parser.add_argument('-N','--fnameTable', help = 'a file of sampl_id and long description for that sample', required=True)
    parser.add_argument('-R','--frange', help = 'a file with the location of exons in the fasta alignment.', required=True)
    parser.add_argument('-A','--f_protein_annotation', help = 'a file with protein_id and annotation. typically a result of protein function annotation', required=True)
    parser.add_argument('-i','--f_ingroupNames', help = 'f_ingroupNames is a file with the target samples. if provided, output alignments with target samples first. default None. if provided, output alignments with target samples first ', default=None)

    f = parser.parse_args()
    prepareRegionsAndAnno(fn=f.fn, faln=f.faln, ftree=f.ftree, fbad=f.fbad, fnameTable=f.fnameTable, frange=f.frange, f_protein_annotation=f.f_protein_annotation, f_ingroupNames=f.f_ingroupNames)
