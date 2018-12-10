# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 17:20:18 2018

@author: ATPs
"""
import numpy as np
from itertools import combinations
from multiprocessing import Pool
import os,sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import mapfileIO
import random
from collections import Counter
import pickle
#A:65 C:67 G:71 T:84 -:45

def countCommon(arr1, arr2):
    '''
    arr1 and arr2 are two numpy array of np.int8
    return the number of common sites
    '''
    return (arr1 == arr2).sum()

def countCommonBases(arr1,arr2):
    '''
    arr1 and arr2 are two numpy array of np.int8
    return a dictionary with the number of common bases and common gaps
    '''
    return {e:((arr1==ord(e)) & (arr2==ord(e))).sum() for e in 'ATCG-'}

def calCommonRatio(arr1, arr2):
    '''
    arr1 and arr2 are two numpy array of np.int8
    return the ratio of common sites
    '''
    return countCommon(arr1,arr2)/len(arr1)

def calCommonNonGapRatio(arr1,arr2):
    '''
    arr1 and arr2 are two numpy array of np.int8, or two filenames
    return the ratio of common non-gap sites
    '''
    if type(arr1) is str:
        arr1 = mapfileIO.loadMapBinary(arr1)
    if type(arr2) is str:
        arr2 = mapfileIO.loadMapBinary(arr2)
    bothNonGap = (arr1 != ord('-')) & (arr2 != ord('-'))
    countBothNonGap = bothNonGap.sum()
    commonBases = countCommon(arr1[bothNonGap], arr2[bothNonGap])
    return commonBases/countBothNonGap

def calCommonNonGapRatio2D(arrs, threads =16):
    '''
    arrs is a list of numpy array of np.int8, or a list of filenames
    return a numpy 2d array with their similarity of matched region
    '''
    arrNum = len(arrs)
    #initiate result
    results = np.zeros([arrNum,arrNum],dtype=float)
    # arr compare with itself, set value 1
    for i in range(arrNum):
        results[i][i] = 1
    
    pairs = list(combinations(range(arrNum),2))
    arr_pairs = [(arrs[e[0]], arrs[e[1]]) for e in pairs]
    pool = Pool(threads)
    pairs_result = pool.starmap(calCommonNonGapRatio,arr_pairs)
    pool.close()
    pool.join()
    for e, v in zip(pairs, pairs_result):
        results[e[0],e[1]] = v
        results[e[1],e[0]] = v
    
    return results

def _choosebase(bases, method=0):
    '''
    helper function of combineNpMap
    choose the best base from a list of int
    '''
    basekeep = [e for e in bases if e != 45]#remove gaps
    if len(basekeep) == 0:# all gap, return gap
        return 45
    if len(set(basekeep)) == 1:#only one kind of base
        return basekeep[0]
    if method == 0:
        return random.choice(basekeep)
    if method == 1:
        return basekeep[0]
    baseCount = Counter(basekeep).most_common()
    if method == 2:
        maxcount = baseCount[0][1]
        basekeep2 = [e[0] for e in baseCount if e[1]==maxcount]
        return random.choice(basekeep2)
    if method == 3:
        mincount = baseCount[-1][1]
        basekeep2 = [e[0] for e in baseCount if e[1]==mincount]
        return random.choice(basekeep2)

def _choosebases(basesNp, method=0):
    '''
    helper function of combineNpMap
    choose the best base from a list of np.array
    '''
    results = []
    for bases in zip(*basesNp):
        results.append(_choosebase(bases,method))
    return np.array(results,dtype=np.int8)


def combineNpMap(files, outfile=None, method=0, threads =16):
    '''
    files is a file stores the location of several binary map files of np.int8
    files can also be a list of file locations, or a list of np.int8 array
    if outfile is None, return np.int8 array
    if outfile is not None, write the result to that file, and return np.int8
    method = 0, random choose each site from the input files
    method = 1, choose bases with the order of input files
    method = 2, choose bases with the highest frequency
    method = 3, choose bases with the lowest frequency
    '''
    #load files
    try:
        arrs = mapfileIO.loadMapBinaries(files,threads=threads)
    except:
        print('input is not file locations, probabally a list of np.array')
        arrs = files
    #split array to threads part
    seqlen = len(arrs[0])
    step = seqlen // threads + 1
    ls_arrs = []
    for start in range(0,seqlen,step):
        arrs_frag = []
        for arr in arrs:
            arrs_frag.append(arr[start:start+step])
        ls_arrs.append(arrs_frag)
    #choose bases
    pool = Pool(threads)
    results = pool.starmap(_choosebases, [[e,method] for e in ls_arrs])
    pool.close()
    pool.join()
    
    result = np.concatenate(results)
    if outfile is not None:
        with open(outfile,'wb') as f:
            pickle.dump(result,f)
            print('done')
    return result
    

def combineMap(files,outfile=None, method=0, threads =16, AsOne=False):
    '''
    files is a file stores the location of several map file 
    files can also be a list with the file locations
    if outfile is None, write to 'output.map' of the working folder
    if outfile is not None, write the result to that file
    only keep 'ATCG'
    method = 0, random choose each site from the input files
    method = 1, choose bases with the order of input files
    method = 2, choose bases with the highest frequency
    method = 3, choose bases with the lowest frequency
    if AsOne is True, the two strand will be used together and output one strand in the output mapfile
    '''
    #prepare outfile
    if outfile is None:
        outfile = 'output.map'
    fout = open(outfile,'w')
    #process the data
    if not AsOne:
        arrs0 = mapfileIO.read2Int8s(files, strand=0, changeGap=True, threads=threads, outfolder=None)
        result0 = combineNpMap(arrs0, outfile=None, method=method, threads =threads)
        del arrs0
        arrs1 = mapfileIO.read2Int8s(files, strand=2, changeGap=True, threads=threads, outfolder=None)
        result1 = combineNpMap(arrs1, outfile=None, method=method, threads =threads)
        del arrs1
        #save the result
        for b0,b1 in zip(result0,result1):
            fout.write(chr(b0)+'\t'+chr(b1)+'\n')
        fout.close()
    else:
        arrs0 = mapfileIO.read2Int8s(files, strand=0, changeGap=True, threads=threads, outfolder=None)
        arrs1 = mapfileIO.read2Int8s(files, strand=2, changeGap=True, threads=threads, outfolder=None)
        arrs = arrs0 + arrs1
        result = combineNpMap(arrs, outfile=None, method=method, threads =threads)
        del arrs0, arrs1, arrs
        #save the result
        fout.write('\n'.join(chr(b) for b in result))
        fout.write('\n')
        fout.close()
    return None

def combineMapSmallMem(filenames, outfile=None, method=0, AsOne=False):
    '''
    work the same as combineMap, single thread, do not load the data into the memory
    '''
    #prepare outfile
    if outfile is None:
        outfile = 'output.map'
    fout = open(outfile,'w')
    
    if type(filenames) is str:
        filenames = open(filenames).readlines()
        filenames = [e.strip() for e in filenames]
    
    ls_files_open = [open(f) for f in filenames]
    for lines in zip(*ls_files_open):
        bases0 = [ord(e[0]) for e in lines if e[0] in 'ATCG']
        bases1 = [ord(e[2]) for e in lines if e[2] in 'ATCG']
        if not AsOne:
            b0 = _choosebase(bases0,method)
            b1 = _choosebase(bases1,method)
            fout.write(chr(b0)+'\t'+chr(b1)+'\n')
        else:
            bases = bases0 + bases1
            b = _choosebase(bases, method)
            fout.write(chr(b)+'\n')
    fout.close()
    print('done')
    return None

    
if __name__ == '__main__':
    pass