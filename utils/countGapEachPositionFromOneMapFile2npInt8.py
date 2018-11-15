
import numpy as np

def countgap(filename,strand=0, outfile = None):
    '''
    given a map file, return a numpy array of 1 or 0 showing whether a position is a gap
    1 means gap, and 0 means 'ATCG'
    strand = 0 or 2, the first or the second strand
    if outfile is None, return the numpy array. Else, save the numpy array in disk
    '''
    if strand == 0:
        strand = 0
    elif strand == 1:
        strand = 2
    else:
        print('wrong strand number!')
        return None
    mapFileGapCout = np.array([e[strand] not in 'ATCG' for e in open(filename)],dtype=np.int8)
    if outfile is None:
        return mapFileGapCout
    #else store the variable in file
    with open(outfile,'wb') as f:
        import pickle
        pickle.dump(mapFileGapCout,f)
        print('done')

description = '''
given a map file, return a numpy array of 1 or 0 showing whether a position is a gap
1 means gap, and 0 means 'ATCG'
strand = 0 or 1, the first or the second strand
if outfile is None, return the numpy array. Else, save the numpy array in disk.
np.int8. Remember to convert to np.int16
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input location of .map files', required=True)
    parser.add_argument('-o','--output',help = 'folder of numpy array of map gap count, default None', required = False, default = None)
    parser.add_argument('-s','--strand', help = 'which strand to use. 0 or 1', default = 0, choices = [0,1], type=int)
    f = parser.parse_args()
    countgap(filename=f.input,strand=f.strand, outfile=f.output)