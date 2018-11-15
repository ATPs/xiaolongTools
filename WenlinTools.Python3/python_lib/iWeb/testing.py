#!/usr/bin/env python

'''
from uniprot import IDmapping

gi='254780193'

ac=['C6XHC3','C6XHC2','C6XHC34']

record=IDmapping(ac,fromdb='ACC+ID',todb='P_GI')

print record.fromid
print record.mapped
print record.toid


from Entrez import IDmapping

gi=['16077368','15966153']

record=IDmapping(gi,fromdb='protein',todb='pubmed',exclusion='protein_pubmed_weighted')

print record.fromid
print record.mapped
print record.toid

'''

from .uniprot import IDsummary

ac=['C6XHC3','C6XHC2','C6XHC34']
record=IDsummary(ac)

print(dir(record))
print(record.esum)


