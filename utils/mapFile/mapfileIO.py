import numpy as np
from multiprocessing import Pool


def fileLine(filename):
    '''
    return line numbers
    '''
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def read2Int8(filename,strand = 0,changeGap=True):
    '''
    filename is a map file, each line with two bases 'A\t\A\n'
    if changGap is True, replace none 'ATCG' to '-'
    strand = 0 or 2, the first or the second strand
    return a numpy array of int8
    '''
    maplines = fileLine(filename)
    mapBaseInt = np.zeros(maplines,dtype=np.uint8) #store bases in np.int8 format
    if strand not in {0,2}:
        print('currently strand can only be 0 or 2')
        return None
    for n,e in enumerate(open(filename)):
        nt = e[strand]
        if changeGap:
            if nt not in 'ATCG':
                nt = '-'
        mapBaseInt[n] = ord(nt)
    return mapBaseInt

def npInt8toSeq(mapBaseInt):
    '''
    bases are coded in a numpy.Int8 array, return a str of bases
    '''
    return ''.join(chr(e) for e in mapBaseInt)

def read2Int8s(filenames, strand=0, changeGap=True, threads=16):
    '''
    filename is list of map files, or a file storing the location of map files
    for each map file, each line with two bases 'A\t\A\n'
    if changGap is True, replace none 'ATCG' to '-'
    strand = 0 or 2, the first or the second strand
    threads is the numbers of CPUs to use
    return a list of numpy array of int8
    '''
    if type(filenames) is str:
        filenames = open(filenames).readlines()
        filenames = [e.strip() for e in filenames]
    pool = Pool(threads)
    mapBaseInts = pool.starmap(read2Int8, [(e, strand, changeGap) for e in filenames])
    pool.close()
    return mapBaseInts

def parse2Int8(filenames, strand = 0, changeGap=True):
    '''
    filename is list of map files, or a file storing the location of map files
    for each map file, each line with two bases 'A\t\A\n'
    if changGap is True, replace none 'ATCG' to '-'
    strand = 0 or 2, the first or the second strand
    return a iterator of numpy array of int8
    '''
    if type(filenames) is str:
        filenames = open(filenames).readlines()
        filenames = [e.strip() for e in filenames]
    for filename in filenames:
        yield read2Int8(filename,strand, changeGap)

def countgap(filename,strand=0, outfile = None):
    '''
    given a map file, return a numpy array of 1 or 0 showing whether a position is a gap, filename can also be a np.Int8 coding the sequences
    1 means gap, and 0 means 'ATCG'
    strand = 0 or 2, the first or the second strand
    if outfile is None, return the numpy array. Else, save the numpy array in disk
    '''
    if strand not in {0,2}:
        print('currently strand can only be 0 or 2')
    if type(filename) is str:
        mapFileGapCout = np.array([e[strand] not in 'ATCG' for e in open(filename)],dtype=np.int8)
    elif type(filename) is np.ndarray:
        mapFileGapCout = np.array([chr(e) not in 'ATCG' for e in filename],dtype=np.int8)
    else:
        print('input filename not np.array, not a valid filename')
        return None
    if outfile is None:
        return mapFileGapCout
    #else store the variable in file
    with open(outfile,'wb') as f:
        import pickle
        pickle.dump(mapFileGapCout,f)
        print('done')