#functions dealing with numpy array object
import numpy as np
import cmn

#output matrix into user-friendly text file
#for example, for 2-d matrix, the result is '1\t2\t0.1234'
def output_matrix(arr,fn=None):
    dimension=arr.ndim
    lines=[]
    for index,content in np.ndenumerate(arr):
        first=''
        for i in range(dimension):
            first+='%s\t' % (index[i])
        lines += ['%s%s' % (first,content)]
    if fn!=None:
        cmn.write_file('\n'.join(lines),fn)
    else:
        return '\n'.join(lines)


def fill_left(arr):
    """
    use the right part to fill in the left part in diagnose matrix
    """
    M,N=arr.shape
    for i in range(M):
        for j in range(i+1,M):#must be diagnosed
            arr[j,i]=arr[i,j]
    return arr


#sum up all the element to make it into one dimension array
#if left empty True, fill in the empty elements
#attetion: when left empty True, matrix must be N*N
def marginal(arr, left_empty=False):
    dim=arr.ndim
    M,N=arr.shape

    try:
        tmp=np.copy(arr)
    except:
        tmp=arr

    if left_empty:
        for i in range(M):
            for j in range(i+1,M):#must be diagnosed
                tmp[j,i]=tmp[i,j]

    for i in range(dim-1):
        tmp=sum(tmp)
    return tmp
