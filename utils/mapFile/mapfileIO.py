import numpy as np
from multiprocessing import Pool
import pickle
import os

def fileLine(filename):
    '''
    return line numbers
    '''
    with open(filename) as f:
        txt = f.read()
        n = txt.count('\n')
        if txt[-1] != '\n':
            n += 1
    return n

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
        if not os.path.exists(outfolder):
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
    
    if isinstance(npFilter,str):
        npFilter = loadMapBinary(npFilter)
    
    if npFilter is None:
        seq = npInt8toSeq(npInt8)
    else:
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
        elif isinstance(npInt8s[0], list) or isinstance(npInt8s[0], tuple):
            print('npInt8s is a list, each element with two elements')
            headers = [e[0] for e in npInt8s]
            npInt8s = [e[1] for e in npInt8s]
        else:
            print('npInt8 is a list, but not in the right format')
            return None
    
    if isinstance(npInt8s,dict):
        print('npInt8 is a dictionary, key is the header for each seq, and value is npInt8 or file of npInt8')
        headers = list(npInt8s.keys())
        npInt8s = list(npInt8s.values())
        
    print('apply headerFun, header', headers[0],'become',headerFun(headers[0]))
    headers = [headerFun(e) for e in headers]
    
    if isinstance(npFilter,str):
        print('npFilter is str, load it from file')
        npFilter = loadMapBinary(npFilter)
    
    if npFilterFun is not None:
        print('try to apply npFilterFun for npFilter')
        if npFilter is None:
            print('npFilter is None, keep all sites')
        else:
            npFilter = npFilterFun(npFilter)
            print('after applying npFilterFun,', npFilter.astype(bool).sum(), 'sites left')
    
    import uuid
    file_shared_variables_in_memory = '/dev/shm/' + str(uuid.uuid4())
    with open(file_shared_variables_in_memory,'wb') as f:
        pickle.dump(npFilter, f)
    npFilter = file_shared_variables_in_memory
    pool = Pool(threads)
    if not merged:
        pool.starmap(npInt8toFasta, [[npInt8, npFilter, header, os.path.join(output,header)] for npInt8,header in zip(npInt8s,headers)])
    if merged:
        fout = open(output,'w')
        parameters = [[npInt8, npFilter, header, None] for npInt8,header in zip(npInt8s,headers)]
        pars = np.array_split(parameters, int(np.ceil(len(parameters)/threads)))
        for par in pars:
            par = [list(e) for e in par]
            seqs = pool.starmap(npInt8toFasta,par)
            for s in seqs:
                fout.write(s)
            seqs = ''
        fout.close()
    pool.close()
    os.remove(npFilter)

def mapFilestoFastas(mapFiles, strand = 0, threads=16, output='fasta', headerFun=lambda x:x.split('_')[0], npFilter=None, npFilterFun=None, merged=True):
    '''
    the same as npInt8stoFastas
    mapFiles is a list of map files
    '''
    print('convert mapFiles to npInt8')
    if type(mapFiles) is str:
        mapFiles = open(mapFiles).readlines()
        mapFiles = [e.strip() for e in mapFiles]
        
    npInt8s = read2Int8s(mapFiles,  strand=strand, changeGap=True, threads=threads, outfolder=None)
    npInt8s = dict(zip(mapFiles, npInt8s))
    print('extract sequence')
    npInt8stoFastas(npInt8s=npInt8s, threads=threads, output=output, headerFun=headerFun, npFilter=npFilter, npFilterFun=npFilterFun, merged=merged)



def npInt8toMultipleFasta(npInt8, header=None, npFilter=None, minlen=0):
    '''
    npInt8 is npInt8 of map file, or a filename of npInt8 file
    npFilter is dictionary with different fragment_id and sites belong to that fragments
    header will be used as the fasta seq name
    return a dictionary with keys in npFilter as key, and a 'fasta' str of nt/AA
    if the str is shorter than minlen, do not include in the final result
    '''
    if header is None:
        header = 'seq'
    
    if isinstance(npFilter,str):
        print('npFilter is str, load it from file')
        npFilter = loadMapBinary(npFilter)
    if isinstance(npFilter,dict):
        print('now npFilter is a dictionary, will output',len(npFilter),'files')
    else:
        print('the format of npFilter is wrong. npFilter type is',type(npFilter))
        return None
    
    if isinstance(npInt8, str):
        npInt8 = loadMapBinary(npInt8)
    if not isinstance(npInt8, np.ndarray):
        print('something wrong with the input format of npInt8, return None')
        return None
    
    dc_seq = {}
    for key in npFilter:
        seq = npInt8toSeq(npInt8[npFilter[key]])
        seqlen = len(seq)
        mlen = minlen
        if mlen <= 1:
            mlen = seqlen * mlen
        if seqlen - seq.count('-') > mlen:
            dc_seq[key] = '>'+header+'\n'+seq+'\n'
        else:
            dc_seq[key] = ''
    
    return dc_seq

def npInt8stoMultipleFastas(npInt8s, threads=16, outputfolder='.', headerFun=lambda x:x.split('_')[0], npFilter=None, minlen=0, minseq = 0):
    '''
    convert npInt8s to fasta sequences. Very similar to npInt8stoFastas. The difference is that here npFilter is a dictionary with key and sites to keep. Thus, the ouput will be multiple files, each with the key as filename
    
    npFilter is a str of binary files of a dictionary, or a dictionary
    
    npInt8s is a filename of file with npInt8 filenames. basename of the npInt8 filenames will be used as header or individual filenames for output
            or a list of npInt8 filenames
            or a list with two elements: header and npInt8 or filename of npInt8
            or a dictionary with header as key, npInt8 or filename of npInt8 as key
    threads is the number of CPUs to run the jobs
    outputfolder is where the files will be stored
    outputfolder should be empty.
    
    minlen is the minimal length of a sequences to be included in the result. Default 0, all sequences will be kept.
    if minlen <= 1: minlen is a ratio of the sequence
    
    minseq is the minimal sequences required to write a file. Default 0, write a file for all fragments in npFilter.
    if minseq <=1: minseq = number_of_sequences * minseq
    '''
    
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
    if callable(headerFun):
        headers = [headerFun(e) for e in headers]

    # make outputfolder if not exist
    if not os.path.exists(outputfolder):
        print('outputfolder does not exist, create one')
        os.makedirs(outputfolder)
    
    #check if outputfolder is empty
    if len(os.listdir(outputfolder)) != 0:
        print('outputfolder is not empty. will not run')
        return None
    
    #create parameters
    ls_params = []
    for npInt8, header in zip(npInt8s, headers):
        ls_params.append([npInt8, header, npFilter, minlen])
    
    pool = Pool(threads)
    ls_results = pool.starmap(npInt8toMultipleFasta, ls_params)
    pool.close()
    
    #minseq filter
    if minseq <= 1:
        minseq = len(npInt8s) * minseq
    #save the results
    dc_results = {key:[] for key in ls_results[0]}
    for result in ls_results:
        for key in result:
            dc_results[key].append(result[key])
    for key, value in dc_results.items():
        if len([e for e in value if e != '']) < minseq:
            continue
        value = ''.join(value)
        if len(value) > 0:
            open(os.path.join(outputfolder,key),'w').write(value)
    
    print('done')

if __name__ == '__main__':
    pass