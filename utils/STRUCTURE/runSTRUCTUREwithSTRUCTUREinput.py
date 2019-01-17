# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 11:11:35 2018

@author: ATPs
"""

import os
import sys
folder_module = os.path.dirname(os.path.realpath(__file__))
#folder_module = r'C:\Users\ATPs\Documents\GitHub\xiaolongTools\utils\STRUCTURE'
sys.path.append(folder_module)
#import fasta2STRUCTUREinput_haploid_nt

default_STRUCTURE = '/home2/s185491/p/STRUCTURE/structure_kernel_src/structure'

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

def runSTRUCTUREwithSTRUCTURE(file_in, NUMINDS, NUMLOCI,SEED='0', MAXPOPS='2', commandsfile = None, STRUCTURE=None, kwargs=''):
    '''
    file_in is a input file to run STRUCTURE directly.
    NUMINDS is the number of individuals in file_in
    NUMLOCI is the number of loci in file_in
    SEED is the seed to run structure, it is a single value or multiple values separated by comma. If it includes multiple values, then multiple structure will run. SEED='0', or SEED='100,298,999'
    MAXPOPS is also a parameter of STRUCTUE which defines the populations for the STRUCTURE run. It can be a single value or or multiple values separated by comma. If it includes multiple values, then multiple structure will run
    kwargs are parameters in STRUCTUE mainparams and extraparams. It should be used as like 'PLOIDY=1,LABEL=1', with no space, key=value, separate by comma
    by default, to be compatible with fasta2STRUCTUREInput_haploid_nt, some default settings include
        PLOIDY = 1
        MAPDISTANCES = 1
        LABEL = 1
        BURNIN = 10000
        NUMREPS = 20000
    for more, check the mainparams and extraparams file in STRUCTURE folder
    STRUCTURE is the location of STRUCTURE program
    only printout the commands to run if commandsfile is None
    '''
    if STRUCTURE is None:
        STRUCTURE = default_STRUCTURE
    
    file_in_folder = os.path.dirname(file_in)
    basename = os.path.basename(file_in)
    #if file_in is inside a folder of "runSTRUCTURE...", use current folder
    if os.path.basename(file_in_folder).startswith('runSTRUCTURE'):
        working_folder = file_in_folder
    else:
        working_folder = os.path.join(file_in_folder,'runSTRUCTURE'+basename)
    if not os.path.exists(working_folder):
        os.makedirs(working_folder)
    
    #change parameters in dc_params based on kwargs
    if len(kwargs) == 0:
        kwargs = {}
    else:
        kwargs = kwargs.split(',')
        kwargs = [e.split('=') for e in kwargs]
        kwargs = {k:v for k,v in kwargs}
    
    for key in kwargs:
        dc_params[key] = kwargs[key]
    dc_params['NUMINDS'] = NUMINDS
    dc_params['NUMLOCI'] = NUMLOCI
    filename_mainparams = os.path.join(working_folder,'mainparams')
    filename_extraparams = os.path.join(working_folder,'extraparams')
    open(filename_mainparams,'w').write(txt_mainparams.format(**dc_params))
    open(filename_extraparams,'w').write(txt_extraparams.format(**dc_params))
    
    #create folders based on SEED and MAXPOPS and cmds
    commandlines = []
    ls_SEED = SEED.split(',')
    ls_MAXPOPS = MAXPOPS.split(',')
    for seed in ls_SEED:
        for maxpop in ls_MAXPOPS:
            outname = basename+'.STRUCTUREresult.s'+seed+'.pop'+maxpop
            commandline = 'cd {working_folder} && {STRUCTURE} -m mainparams -e extraparams -K {MAXPOPS} -L {NUMLOCI} -N {NUMINDS} -i {file_in} -o {file_out} -D {SEED}'.format(working_folder=working_folder, STRUCTURE=STRUCTURE, MAXPOPS=maxpop, NUMLOCI=NUMLOCI, NUMINDS=NUMINDS, file_in=file_in, file_out = outname, SEED=seed)
            #print(commandline)
            commandlines.append(commandline)
        
    if commandsfile is not None:
        with open(commandsfile,'a') as f:
            for commandline in commandlines:
                f.write(commandline+'\n')
    
    for commandline in commandlines:
        print(commandline)
    print('\ntotally', len(commandlines), 'commands to run')
    
    return commandlines


description = '''file_in is a input file to run STRUCTURE directly.
    NUMINDS is the number of individuals in file_in
    NUMLOCI is the number of loci in file_in
    SEED is the seed to run structure, it is a single value or multiple values separated by comma. If it includes multiple values, then multiple structure will run. SEED='0', or SEED='100,298,999'
    MAXPOPS is also a parameter of STRUCTUE which defines the populations for the STRUCTURE run. It can be a single value or or multiple values separated by comma. If it includes multiple values, then multiple structure will run
    kwargs are parameters in STRUCTUE mainparams and extraparams. It should be used as like 'PLOIDY=1,LABEL=1', with no space, key=value, separate by comma
    by default, to be compatible with fasta2STRUCTUREInput_haploid_nt, some default settings include
        PLOIDY = 1
        MAPDISTANCES = 1
        LABEL = 1
        BURNIN = 10000
        NUMREPS = 20000
    for more, check the mainparams and extraparams file in STRUCTURE folder
    STRUCTURE is the location of STRUCTURE program
    only printout the commands to run if commandsfile is None
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--file_in', help = 'a input file to run STRUCTURE directly.', required=True)
    parser.add_argument('-I','--NUMINDS',help = 'the number of individuals in file_in', required=True)
    parser.add_argument('-L','--NUMLOCI',help = 'the number of loci in file_in', required=True)
    parser.add_argument('-s','--SEED',help = '''the seed to run structure, it is a single value or multiple values separated by comma. If it includes multiple values, then multiple structure will run. SEED='0', or SEED='100,298,999' default=0''', default = 0)
    parser.add_argument('-P','--MAXPOPS',help = 'a parameter of STRUCTUE which defines the populations for the STRUCTURE run. It can be a single value or or multiple values separated by comma. If it includes multiple values, then multiple structure will run. default=2',default = 2)
    parser.add_argument('-o','--commandsfile',help = 'save commandlines to this file. default None, only print the commandlines',default = None)
    parser.add_argument('-S','--STRUCTURE',help = 'location of STRUCTURE program, default None, use the one in bioHPC',default = None)
    parser.add_argument('-k','--kwargs',help = '''kwargs are parameters in STRUCTUE mainparams and extraparams. It should be used as like 'PLOIDY=1,LABEL=1', with no space, key=value, separate by comma''',default = '')
    f = parser.parse_args()
    runSTRUCTUREwithSTRUCTURE(file_in=f.file_in, NUMINDS=f.NUMINDS, NUMLOCI=f.NUMLOCI,SEED=f.SEED, MAXPOPS=f.MAXPOPS, commandsfile = f.commandsfile, STRUCTURE=f.STRUCTURE, kwargs=f.kwargs)