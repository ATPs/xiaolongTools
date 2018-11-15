

from multiprocessing import Pool
import numpy as np
import os
import pickle

file_shared_variables_in_memory = '/dev/shm/shared_variables_in_memory'

def filename2seqname(filename):
    '''
    filename of map file is like '5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map' or '/archive/butterfly/maps/debiased/5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map'
    return '5518' as the sequence name
    '''
    return os.path.basename(filename).split('_')[0]

def extractSequenceFromMapFileInSites(filename,strand = 0):
    '''
    filename is a map file, each line with two bases 'A\t\A\n'
    replace none 'ATCG' to '-'
    sites is list, first element is fragment_ID, second element is a set of int, with sites in the map file that will be kept.
    strand = 0 or 2, the first or the second strand
    return a list, each element with fragment_ID, sequence
    fragment_ID will be used as output filename
    sequence is in fasta format. seqname is generate by filename2seqname
    '''
    dc_nt = {}
    with open(file_shared_variables_in_memory,'rb') as f:
        sites, sites_all = pickle.load(f)
    #sites_all = set.union(*[e[1] for e in sites])
    for n,e in enumerate(open(filename)):
        n = str(n)
        nt = e[strand]
        if nt not in 'ATCG':
            nt = '-'
        if n in sites_all:
            dc_nt[n] = nt
    
    #print(filename,'dc_nt len',len(dc_nt))
    del sites_all
    seqname = filename2seqname(filename=filename)
    fragmentID_sequence = []
    for fragment_ID, fragment_sites in sites:
        seq = ''.join([dc_nt[e] for e in fragment_sites])
        if len(seq) != seq.count('-'):
            fragmentID_sequence.append([fragment_ID,'>'+seqname+'\n'+seq+'\n'])
    
    return fragmentID_sequence

def extractSequencesFromMapFilesInSites(filenames,strand = 0, threads = 16):
    '''
    given a list of map files, each line with two bases 'A\t\A\n'
    sites is a set of int, with sites in the map file that will be kept.
    strand = 0 or 2, the first or the second strand
    threads, number of thread to use
    return a list of sequence of nucleotides
    '''
    #sites_all = set.union(*[e[1] for e in sites])
    pool = Pool(threads)
    fragmentID_sequences = pool.starmap(extractSequenceFromMapFileInSites, [(e,strand) for e in filenames])
    pool.close()
    
    #store the result in a dictionary. key is fragment_ID, value is fasta sequences in list
    dc_fragmentID2sequences = {}
    for fragmentID_sequence in fragmentID_sequences:
        for fragment_ID, sequence in fragmentID_sequence:
            if fragment_ID not in dc_fragmentID2sequences:
                dc_fragmentID2sequences[fragment_ID] = []
            dc_fragmentID2sequences[fragment_ID].append(sequence)
    
    #save files
    for fragment_ID, sequences in dc_fragmentID2sequences.items():
        open(fragment_ID,'a').write(''.join(sequences))


def mapfilelineprocessing(x):
    x1,x2 = x.split(':')
    x2 = tuple(x2.split()) # use tuple to save memory
    return x1,x2

def getSequencesFromMapFilesInSites(files_map,file_sites,outfolder,strand = 0, threads = 16):
    '''
    files_map is a file store the location of map files, or a, or a list of map filepath
    map file, each line with two bases 'A\t\A\n'
    mapfilename is like '5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map', sequence name will be '5518'
    file_sites is a file, or a set of int, with sites in the map file that will be kept.
    outfolder is where the output files will be stored.
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
    
    if type(file_sites) is list:
        sites = file_sites
    else:
        sites_lis = open(file_sites).readlines()
        sites = list(map(mapfilelineprocessing,sites_lis))
    print('number of output files will be',len(sites))
    
    # process the files in batch
    if os.path.exists(outfolder):
        import shutil
        print('outfolder exist, everything in it will be removed!')
        shutil.rmtree(outfolder)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    os.chdir(outfolder)
    
    sites_all = set.union(*[set(e[1]) for e in sites])
    with open(file_shared_variables_in_memory,'wb') as f:
        pickle.dump([sites, sites_all],f)
        del sites
        del sites_all
    file_batches = [files[i:min(i+threads, len(files))] for i in range(0,len(files),threads)]
    print('total batches',len(file_batches))
    threads = int(len(files)/len(file_batches)) + 1
    file_batches = [files[i:min(i+threads, len(files))] for i in range(0,len(files),threads)]
    print('re-adjust thread to balance the work load. set threads to ', threads)
    for n, files_processing in enumerate(file_batches):
        print('processing batch',n+1)
        extractSequencesFromMapFilesInSites(filenames=files_processing, strand = strand, threads = threads)
    
    print('finished!')
    if os.path.exists(file_shared_variables_in_memory):
        os.remove(file_shared_variables_in_memory)



description = '''input a file with filepath map files and sites to keep for different fragments, output multiple fasta files to the output folder. This program can also be imported in other modules to return the value directly.

The sites to keep in stored in the format
        fragment_ID1:123 124 ... 125(line number in the .map file)
        fragment_ID2:123 124 ... 125(line number in the .map file)

output multiple files to the output folder. sequence names is processed as:
    filename of map file is like '5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map' or '/archive/butterfly/maps/debiased/5518_Junonia_coenia_JC_v1.0.scaffolds_snp_step2.map'
    return '5518' as the sequence name

filename will be fragment_IDs

Note: empty sequences will be removed (sequences with only '-')
      outfolder should be empty/new, otherwise everything in it will be removed.
      be careful about the memory consumption. If the input uses a lot of memory, increase the thread number won't help much. 16 is a good number.
'''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file storing the location of .map files', required=True)
    parser.add_argument('-t','--threads',help = 'number of threads to use', default = 16,type=int)
    parser.add_argument('-b','--sites',help = 'sites that will be extracted. ', required = True)
    parser.add_argument('-o','--output',help = 'folder of the where the output files stored', required = True)
    parser.add_argument('-s','--strand', help = 'which strand to use. 0 or 1', default = 0, choices = [0,1], type=int)
    f = parser.parse_args()
    getSequencesFromMapFilesInSites(files_map=f.input, file_sites=f.sites, strand=f.strand, outfolder=f.output, threads=f.threads)
    

