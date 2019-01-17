# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 12:30:47 2018

@author: ATPs
"""

import sys
filein = sys.argv[1]

lines = open(filein).readlines()
for line in lines:
    line=line.strip()
    words = line.split('\t')
    if len(words) <8:
        print(line)
        continue
    anno = words[8]
    anno_ls = anno.split(';')
    anno_ls = [e.strip().split() for e in anno_ls]
    if words[2] == 'transcript':
        words[2] = 'mRNA'
        anno_new = 'ID='+anno_ls[1][1].replace('"','')#+';geneID='+anno_ls[0][1].replace('"','')
    if words[2] == 'exon':
        words[2] = 'cds'
        anno_new = 'Parent='+anno_ls[1][1].replace('"','')#+';geneID='+anno_ls[0][1].replace('"','')
    words[8] = anno_new
    print('\t'.join(words))
    