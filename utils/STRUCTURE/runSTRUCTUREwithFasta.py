# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 11:11:35 2018

@author: ATPs
"""

import os
import sys
folder_module = os.path.realpath(__file__)
#folder_module = r'C:\Users\ATPs\Documents\GitHub\xiaolongTools\utils\STRUCTURE'
sys.path.append(folder_module)
import fasta2STRUCTUREinput_haploid_nt

STRUCTURE = '/home2/s185491/p/STRUCTURE/structure_kernel_src/structure'

file_mainparams = os.path.join(folder_module,'mainparams.ref')
file_mainparams_default = os.path.join(folder_module,'mainparams')
file_extraparams = os.path.join(folder_module,'extraparams.ref')
file_extraparams_default = os.path.join(folder_module,'extraparams')

txt_mainparams = open(file_mainparams).read()
txt_extraparams = open(file_extraparams).read()

def get_params_default(file_default):
    '''
    return a dictionary with the parameters that need to be determined in file_params
    with the default value
    '''
    ls_txt = open(file_default).readlines()
    ls_txt = [e for e in ls_txt if e.startswith('#define')]
    dc = {}
    for l in ls_txt:
        es = l.split()
        dc[es[1]] = es[2]
    return dc

dc_mainparams = get_params_default(file_mainparams_default)
dc_extraparams = get_params_default(file_extraparams_default)
dc_params = dc_mainparams.copy()
dc_params.update(dc_extraparams)

def runSTRUCTUREwithSTRUCTURE(file_in, NUMINDS, NUMLOCI,SEED='0', MAXPOPS='2', STRUCTURE=STRUCTURE, **kwargs):
    '''
    file_in is a input file to run STRUCTURE directly.
    NUMINDS is the number of individuals in file_in
    NUMLOCI is the number of loci in file_in
    SEED is the seed to run structure, it is a single value or multiple values separated by comma. If it includes multiple values, then multiple structure will run
    MAXPOPS is also a parameter of STRUCTUE which defines the populations for the STRUCTURE run. It can be a single value or or multiple values separated by comma. If it includes multiple values, then multiple structure will run
    '''
    file_in_folder = os.path.dirname(file_in)
    basename = os.path.basename(file_in)
    #if file_in is inside a folder of "runSTRUCTURE...", use current folder
    if os.path.basename(file_in_folder).startswith('runSTRUCTURE'):
        working_folder = file_in_folder
    else:
        working_folder = os.path.join(file_in_folder,'runSTRUCTURE'+basename)
    if not os.path.exists(working_folder):
        os.makedirs(working_folder)

def runSTRUCTURE(file_in, isfasta=True, SEED='0', MAXPOPS='2', **kwargs):
    '''
    file_in is a input file to run STRUCTURE. 
    If isfasta=True, convert it to the input of STRUCTURE and run structures
    else, file_in is a input for STRUCTURE, no conversion. But the NUMINDS and NUMLOCI need to be provided in kwargs
    SEED is the seed to run structure, it is a single value or multiple values separated by comma. If it includes multiple values, then multiple structure will run
    MAXPOPS is also a parameter of STRUCTUE which defines the populations for the STRUCTURE run. It can be a single value or or multiple values separated by comma. If it includes multiple values, then multiple structure will run
    **kwargs are parameters in STRUCTUE mainparams and extraparams
    '''
    file_in_folder = os.path.dirname(file_in)
    basename = os.path.basename(file_in)
    #if file_in is inside a folder of "runSTRUCTURE...", use current folder
    if os.path.basename(file_in_folder).startswith('runSTRUCTURE'):
        working_folder = file_in_folder
    else:
        working_folder = os.path.join(file_in_folder,'runSTRUCTURE'+basename)
    if not os.path.exists(working_folder):
        os.makedirs(working_folder)
    
    #prepare the input and output filenames
    INFILE = basename
    seeds = SEED.split(',')
    maxpops = MAXPOPS.split(',')
    for s in seeds:
        for m in maxpops:
            
    
    if isfasta:
        print('input file_in is a fasta file, convert it to the input of STRUCTURE')
        fasta2STRUCTUREinput_haploid_nt(filename,outfile = 'default', gapcut = 0.8, Ncut=4, missing_data = -9, threads = 32, header_map_dist=True)
        
        
        