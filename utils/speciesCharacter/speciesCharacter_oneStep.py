# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 16:27:41 2019

@author: ATPs
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import compute_conserved_letter_dna_withCloseGroup_exclusion_threads
import prepare_regions_and_anno

description = """
two functions in one step
def computeConservedLetterDNAwithCloseGroupExclusion(fa, fn, fclose, fbad, fout, fexclude, frange, frename, gapInGroup=1.0, gapInClose=0.9, gapOutGroup=0.7, consOutGroup=0.7, outfilename = None):
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


def prepareRegionsAndAnno(fn, faln, ftree, fbad, fnameTable, frange, f_protein_annotation):
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
    '''
    
    
"""


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
    parser.add_argument('-t','--ftree', help = 'a tree file which will be used to re-order samples')
    parser.add_argument('-A','--f_protein_annotation', help = 'a file with protein_id and annotation. typically a result of protein function annotation')
    parser.add_argument('-T','--threads', help = '''number of CPUs to use. default 16''', default=16, type=int)
    f = parser.parse_args()
    print('process step1')
    compute_conserved_letter_dna_withCloseGroup_exclusion_threads.computeConservedLetterDNAwithCloseGroupExclusion(fa=f.fa, fn=f.fn, fclose=f.fclose, fbad=f.fbad, fout=f.fout, fexclude=f.fexclude, frange=f.frange, frename=f.frename, gapInGroup=f.gapInGroup, gapInClose=f.gapInClose, gapOutGroup=f.gapOutGroup, consOutGroup=f.consOutGroup, outfilename=f.outfilename, threads=f.threads)
    print('step1 finished')
    if f.outfilename is None:
        filename = '%s_%s_correctIndex.cc' % (os.path.basename(f.fa), os.path.basename(f.fn))
    else:
        filename = f.outfilename
    prepare_regions_and_anno.prepareRegionsAndAnno(fn=filename, faln=f.fa, ftree=f.ftree, fbad=f.fbad, fnameTable=f.frename, frange=f.frange, f_protein_annotation=f.f_protein_annotation,f_ingroupNames=f.fn)
