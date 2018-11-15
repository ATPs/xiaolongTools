#!/bin/env python

#change the parameter style from old to new

#author:wenlin li; date: 2011-10-13

import re
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
def reformat(command):
    m=re.compile('.*\.py ')#find the filename

    #erase filename
    k=m.search(command)
    word=k.group()
    command=command.replace(word,'')

    #-I T
    if '-I T' not in command:
        command+=' -show_gis'
    else:
        command=command.replace('-I T','-show_gis')

    #program name
    blastcmd = "/home2/wli/local/ncbi-blast-2.2.29+/bin/psiblast "
    if 'blastpgp' in command:
        command=command.replace('blastpgp',blastcmd)
    else:
        command=blastcmd + command
    #-b and -v
    if '-b ' not in command:
        command+=' -b 5000 '
    if '-v ' not in command:
        command+=' -v 5000 '

    #-d
    if '-d' not in command:
        command+=' -db /usr2/db/fasta/nr '
    else:
        command=command.replace('-d','-db')

    #-m
    if '-m' not in command:
        command+=' -m 0 '

    command=command.replace('-cpu','-num_threads')
    command=command.replace('-i','-query')
    command=command.replace('-o','-out')
    command=command.replace('-j ','-num_iterations ')
    command=command.replace('-b ','-num_alignments ')
    command=command.replace('-v ','-num_descriptions ')
    command=command.replace('-C ','-out_pssm ')
    command=command.replace('-R ','-in_pssm ')
    command=command.replace('-B ','-in_msa ')
    command=command.replace('-e ','-evalue ')
    command=command.replace('-h ','-inclusion_ethresh ')
    command=command.replace('-m ','-outfmt ')
    return command

if __name__=='__main__':
    import sys,os
    if len(sys.argv)==1:#no parameter
        print('use this program as blastpgp!\n')
        sys.exit()

    #read record
    command=''
    for i in sys.argv:
        command+=i+' '
    cmd=reformat(command)
    print('running command: '+cmd)
    os.system(cmd)

