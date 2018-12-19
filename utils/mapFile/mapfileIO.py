import numpy as np
from multiprocessing import Pool
import pickle
import os

def fileLine(filename):
    '''
    return line numbers
    '''
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1

def read2Int8(filename,strand = 0,changeGap=True, outfile=None):
    '''
    filename is a map file, each line with two bases 'A\t\A\n'
    if changGap is True, replace none 'ATCG' to '-'
    strand = 0 or 2, the first or the second strand
    if outfile is None, do not output a file, else, save map file in binary format file
    else, return a numpy array of int8
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
    if outfile is not None:
        with open(outfile,'wb') as f:
            pickle.dump(mapBaseInt,f)
            print('done')
    else:
        return mapBaseInt

def npInt8toSeq(mapBaseInt):
    '''
    bases are coded in a numpy.Int8 array, return a str of bases
    '''
    return ''.join(chr(e) for e in mapBaseInt)


def read2Int8s(filenames, strand=0, changeGap=True, threads=16, outfolder=None):
    '''
    filename is list of map files, or a file storing the location of map files
    for each map file, each line with two bases 'A\t\A\n'
    if changGap is True, replace none 'ATCG' to '-'
    strand = 0 or 2, the first or the second strand
    threads is the numbers of CPUs to use
    return a list of numpy array of int8
    if outfolder is not None, save map file in binary format file, filename is the same as input files
    '''
    if type(filenames) is str:
        filenames = open(filenames).readlines()
        filenames = [e.strip() for e in filenames]
    
    pool = Pool(threads)
    if outfolder is not None:
        os.makedirs(outfolder)
    
    if outfolder is None:
        mapBaseInts = pool.starmap(read2Int8, [(e, strand, changeGap,None) for e in filenames])
        pool.close()
        pool.join()
        return mapBaseInts
    else:
        pool.starmap(read2Int8, [(e, strand, changeGap, os.path.join(outfolder,os.path.basename(e))) for e in filenames])
        pool.close()
        pool.join()
        return None

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
        pickle.dump(mapFileGapCout,f)
        print('done')

def loadMapBinary(filename):
    '''
    filename is the binary format of np.int8 array, load it from the disk
    '''
    with open(filename,'rb') as f:
        return pickle.load(f)

def loadMapBinaries(filenames, threads=16):
    '''
    filename is list of map files, or a file storing the location of map files
    return a list of np array with the sequences coding in np.int8
    '''
    if type(filenames) is str:
        filenames = open(filenames).readlines()
        filenames = [e.strip() for e in filenames]
    if type(filenames) is not list:
        print('wrong input format')
        return None
    pool = Pool(threads)
    results = pool.map(loadMapBinary, filenames)
    pool.close()
    pool.join()
    return results

def npInt8toFasta(npInt8, npFilter=None, header=None,outfile=None):
    '''
    npInt8 is npInt8 of map file, or a filename of npInt8 file
    npFilter is numpy array the same length with npInt8 or a file of numpy array binary which determines sites to keep in npInt8
        if npFilter is None, keep all sites
    if outfile is None, return a str of fasta format
    else, write fasta seq to outfile
    if header is None:
        if outfile is not None:
            header = basename of outfile
        else:
            header = 'seq'
    '''
    if header is None:
        if outfile is not None:
            header = os.path.basename(outfile)
        else:
            header = 'seq'
    
    if isinstance(npInt8, str):
        npInt8 = loadMapBinary(npInt8)
    
    if not isinstance(npInt8, np.ndarray):
        print('something wrong with the input format of npInt8, return None')
        return None
    
    if npFilter is None:
        seq = npInt8toSeq(npInt8)
    else:
        if isinstance(npFilter,str):
            npFilter = loadMapBinary(npFilter)
        if not isinstance(npFilter, np.ndarray):
            print('something wrong with the input format of npFilter, return None')
            return None
        seq = npInt8toSeq(npInt8[npFilter])
    
    fasta = '>'+header+'\n' + seq+'\n'
    if outfile is None:
        return fasta
    with open(outfile,'w') as f:
        f.write(fasta)
    return None

def npInt8stoFastas(npInt8s, threads=16, output='fasta', headerFun=lambda x:x.split('_')[0], npFilter=None, npFilterFun=None, merged=True):
    '''
    npInt8s is a filename of file with npInt8 filenames. basename of the npInt8 filenames will be used as header or individual filenames for output
            or a list of npInt8 filenames
            or a list with two elements: header and npInt8 or filename of npInt8
            or a dictionary with header as key, npInt8 or filename of npInt8 as key
    threads is the number of CPUs to run the jobs
    output is where to store the output. 
        if merged is True, out is the filename to store the output
        else, out is the outfolder to store the output files
    headerFun is the function applied to headers for each sequence, default keep the part before '_'
    npFilter is the numpy array to filter the input npInt8s, it can a str or np.array
    npFilterFun is the function that will be applied to npFilter, to further process the fulter parameter
    merged is True, output the sequences to a single file, else, to multiple files
    '''
    if merged:
        print('output will be merged to a single file')
    else:
        print('output each sequence to an individual file')
    
    if isinstance(npInt8s, str):
        print('npInt8s is a filename of file with npInt8 files')
        npInt8s = open(npInt8s).readlines()
        npInt8s = [e.strip() for e in npInt8s]
        headers = [os.path.basename(e) for e in npInt8s]
    
    if isinstance(npInt8s, list):
        if isinstance(npInt8s[0], str):
            print('npInt8s is a list of npInt8 filenames')
            headers = [os.path.basename(e) for e in npInt8s]
        elif isinstance(npInt8s[0], list):
            print('npInt8s is a list, each element with two elements')
            headers = [e[0] for e in npInt8s]
            npInt8s = [e[1] for e in npInt8s]
        else:
            print('npInt8 is a list, but not in the right format')
    
    if isinstance(npInt8s,dict):
        print('npInt8 is a dictionary, key is the header for each seq, and value is npInt8 or file of npInt8')
        headers = list(npInt8s.keys())
        npInt8s = list(npInt8s.values())
        
    print('apply headerFun, header', headers[0],'become',headerFun(headers[0]))
    headers = [headerFun(e) for e in headers]
    
    if npFilterFun is not None:
        print('try to apply npFilterFun for npFilter')
        if npFilter is None:
            print('npFilter is None, keep all sites')
        else:
            npFilter = npFilterFun(npFilter)
            print('after applying npFilterFun,', npFilter.astype(bool).sum(), 'sites left')
    
    pool = Pool(threads)
    if not merged:
        pool.starmap(npInt8toFasta, [[npInt8, npFilter, header, os.path.join(output,header)] for npInt8,header in zip(npInt8s,headers)])
    if merged:
        fout = open(output,'w')
        parameters = [[npInt8, npFilter, header, None] for npInt8,header in zip(npInt8s,headers)]
        pars = np.array_split(parameters, int(np.ceil(len(parameters)/threads)))
        for par in pars:
            seqs = pool.starmap(npInt8toFasta,par)
            for s in seqs:
                fout.write(s)
            seqs = ''
        fout.close()
    pool.close()
    

if __name__ == '__main__':
    pass