

from multiprocessing import Pool
import numpy as np

def countgap(filename,strand=0):
    '''
    given a map file, return a numpy array of 1 or 0 showing whether a position is a gap
    1 means gap, and 0 means 'ATCG'
    strand = 0 or 2, the first or the second strand
    '''
    return np.array([e[strand] not in 'ATCG' for e in open(filename)],dtype=np.int16)

def countgaps(filenames,strand = 0, threads = 16):
    '''
    given a list of map filenames, a numpy array the number of gaps in each position.
    those map filenames should be from the same genome
    strand = 0 or 3, the first or the second strand
    '''
    pool = Pool(threads)
    gaps = pool.map(countgap,filenames)
    pool.close()
    return sum(gaps)

def array2str(np_array):
    '''
    convert a numpy array to str, joined by '\n'
    '''
    np_array = np_array.astype(str)
    return '\n'.join(np_array)

def array2str_thread(np_array,threads = 16):
    '''
    convert a numpy str array to str, joined by '\n'
    with thread enabled
    '''
    pool = Pool(threads)
    substrs = pool.map(array2str, np.array_split(np_array,threads))
    pool.close()
    return '\n'.join(substrs)

def gapCountsConvert2str(gapcounts,threads = 16):
    '''
    gapcounts is the np.int16 array. 
    convert it to a long str joined by '\n' to write with threads to speed up
    '''
#    gapstr = gapcounts.astype(str)
    gapstr = gapcounts
    gapstrs = np.array_split(gapstr,10)
    substrs = []
    for s in gapstrs:
        substrs.append(array2str_thread(s,threads=threads))
    return '\n'.join(substrs)
    

def countGapEachPositionFromMapFiles(files_map,strand=0,outfile = None,threads = 16):
    '''
    files_map is a file of a list of .map files, each with lines of nucleotide in two strands
    or files_map is a python list of .map files
    threads is the number of threads to use. Default 16
    strand = 0 or 1, the first or the second strand
    return a numpy array of np.int16 with counts of gaps in each position
    if outfile is not None, write a text file to store the result
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
    
    # to save memories, apply the function in batch
    file_batches = [files[i:min(i+threads, len(files))] for i in range(0,len(files),threads)]
    gapcounts = countgaps(file_batches[0], strand=strand,threads=threads)
    l_bases = len(gapcounts)
    print('length of map file is', l_bases)
    for n in range(1,len(file_batches)):
        print('procession',threads * n)
        _gaps = countgaps(file_batches[n], strand=strand,threads=threads)
        gapcounts += _gaps
    
    # save the file
    if outfile is not None:
        with open(outfile,'w') as fout:
            fout.write(gapCountsConvert2str(gapcounts=gapcounts,threads=threads))
            fout.write('\n')
    
    return gapcounts

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input file of location of many map files, output file with the count of gaps in each base location. This program can also be imported in other modules to return the value directly')
    parser.add_argument('-i','--input', help = 'input file storing the location of .map files', required=True)
    parser.add_argument('-t','--threads',help = 'number of threads to use', default = 16,type=int)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored', required = True)
    parser.add_argument('-s','--strand', help = 'which strand to use. 0 or 1', default = 0, choices = [0,1], type=int)
    f = parser.parse_args()
    countGapEachPositionFromMapFiles(files_map=f.input, strand=f.strand, outfile=f.output, threads=f.threads)

