

from multiprocessing import Pool
import numpy as np
import os
import pickle

#store the passed sets in memory, which can be pretty large to save time
file_shared_variables_in_memory = '/dev/shm/shared_variables_in_memory_sites_of_set_temp'

def extractSequenceFromMapFileInSites(filename,strand = 0):
    '''
    filename is a map file, each line with two bases 'A\t\A\n'
    replace none 'ATCG' to '-'
    sites is a set of int, with sites in the map file that will be kept.
    strand = 0 or 2, the first or the second strand
    return a sequence of nucleotides
    '''
    s = []
    with open(file_shared_variables_in_memory,'rb') as f:
        sites = pickle.load(f)
    for n,e in enumerate(open(filename)):
        if n in sites:
            nt = e[strand]
            if nt not in 'ATCG':
                nt = '-'
            s.append(nt)
    return ''.join(s)

def extractSequencesFromMapFilesInSites(filenames,strand = 0, threads = 16):
    '''
    given a list of map files, each line with two bases 'A\t\A\n'
    sites is a set of int, with sites in the map file that will be kept.
    strand = 0 or 2, the first or the second strand
    threads, number of thread to use
    return a list of sequence of nucleotides
    '''
    pool = Pool(threads)
    seqs = pool.starmap(extractSequenceFromMapFileInSites, [(e,strand) for e in filenames])
    pool.close()
    return seqs

def filename2seqname(filename):
    '''
    filename of map file is like '5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map' or '/archive/butterfly/maps/debiased/5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map'
    return '5518' as the sequence name
    '''
    return os.path.basename(filename).split('_')[0]

def getSequencesFromMapFilesInSites(files_map,file_sites,outfile=None,strand = 0, threads = 16):
    '''
    files_map is a file store the location of map files, or a, or a list of map filepath
    map file, each line with two bases 'A\t\A\n'
    mapfilename is like '5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map', sequence name will be '5518'
    file_sites is a file, or a set of int, with sites in the map file that will be kept.
    if outfile is None, do not output fasta file, return a txt file in fasta format.
    if outfile is not None, write the fasta file to outfile.
    strand = 0 or 1, the first or the second strand
    threads, number of thread to use
    '''
    if strand == 0:
        strand = 0
    elif strand == 1:
        strand = 2
    else:
        print('wrong strand number!')
        return None
    print('threads used: ', threads)
    
    if type(files_map) is list:
        files = files_map
    else:
        files = open(files_map).read().split('\n')
        files = [e for e in files if e != ''] # remove empty files
    print('number of map files is ', len(files))
    
    if type(file_sites) is set:
        sites = file_sites
    else:
        sites = open(file_sites).read().split()
        sites = set([int(e) for e in sites])
    print('number of sites to keep is',len(sites))
    with open(file_shared_variables_in_memory,'wb') as f:
        pickle.dump(sites,f)
        del sites
    
    
    # process the files in batch
    if outfile is None:
        l_seqs = []
    else:
        fout = open(outfile,'w')
    file_batches = [files[i:min(i+threads, len(files))] for i in range(0,len(files),threads)]
    print('total batches',len(file_batches))
    threads = int(len(files)/len(file_batches)) + 1
    file_batches = [files[i:min(i+threads, len(files))] for i in range(0,len(files),threads)]
    print('re-adjust thread to balance the work load. set threads to ', threads)
    
    for n, files_processing in enumerate(file_batches):
        print('processing batch',n+1)
        seqs = extractSequencesFromMapFilesInSites(filenames=files_processing, strand = strand, threads = threads)
        seqs_names = [filename2seqname(f) for f in files_processing]
        seqs_fasta = ['>'+seqs_names[n]+'\n'+seqs[n]+'\n' for n in range(len(files_processing))]
        if outfile is None:
            l_seqs = l_seqs + seqs_fasta
        else:
            fout.write(''.join(seqs_fasta))
            del seqs
            del seqs_fasta
    
    if outfile is None:
        return l_seqs
    else:
        fout.close()
    
    print('finished!')
    if os.path.exists(file_shared_variables_in_memory):
        os.remove(file_shared_variables_in_memory)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a file with filepath map files and sites to keep, output a fasta file. This program can also be imported in other modules to return the value directly')
    parser.add_argument('-i','--input', help = 'input file storing the location of .map files', required=True)
    parser.add_argument('-t','--threads',help = 'number of threads to use. default 16', default = 16,type=int)
    parser.add_argument('-b','--sites',help = 'sites that will be extracted', required = True)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored', required = True)
    parser.add_argument('-s','--strand', help = 'which strand to use. 0 or 1. default 0', default = 0, choices = [0,1], type=int)
    f = parser.parse_args()
    getSequencesFromMapFilesInSites(files_map=f.input, file_sites=f.sites, strand=f.strand, outfile=f.output, threads=f.threads)

