#!/usr/bin/env python3

from bs4 import BeautifulSoup
import requests
gene_id=1

def getPubmedIDs(gene_id):
    '''
    with a gene_id like 1
    return a list of pubmed_ids
    '''
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={gene_id}&retmode=xml'
    r = requests.get(url)
    soup = BeautifulSoup(r.text)
    pmids = soup.findAll('pubmedid')
    pmids = [e.text for e in pmids]
    return len(set(pmids))

import pandas as pd
import time
df = pd.read_csv('/lab02/Data_Raw/Xiaolong/HumanCommon/HGNC/hgnc_complete20190821.txt',sep='\t')
df = df[['HGNC ID', 'Approved symbol', 'NCBI Gene ID', 'Ensembl gene ID']]
df = df[df['NCBI Gene ID'].notna()]
df = df.sort_values(by='NCBI Gene ID')

for r, row in df.iterrows():
    gene_id = int(row['NCBI Gene ID'])
    try:
        pmid_count = getPubmedIDs(gene_id)
        df.loc[r,'pmid_count'] = pmid_count
        print(gene_id, pmid_count)
    except:
        print(gene_id, 'not finished')
    time.sleep(0.125)

df['NCBI Gene ID'] = df['NCBI Gene ID'].astype(int)
df['pmid_count'] = df['pmid_count'].fillna(0)
df['pmid_count'] = df['pmid_count'].astype(int)
df.to_csv('/lab02/Data_Raw/Xiaolong/HumanCommon/HGNC/hgnc_complete20190821gene_pmidCount.txt',sep='\t', index=None)

df2 = df.sort_values(by='pmid_count',ascending=False)
df2 = df2.iloc[:200]
df2.to_csv('/lab02/Data_Raw/Xiaolong/HumanCommon/HGNC/hgnc_complete20190821gene_pmidCountTop200.txt',sep='\t', index=None)


