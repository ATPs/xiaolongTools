
from ete3 import Tree

def orderTree(tree, outgroup = 'J_neildi'):
    '''
    tree is a text of tree in newick format or a ete3 tree
    return a text of tree in newick format
    set outgroup, and order tree nodes
    '''
    t = Tree(';')
    if type(tree) is str:
        tree = Tree(tree)
    elif type(tree) == type(t):
        tree = tree
    else:
        print('wrong format', tree)
        return None
    tree.set_outgroup(outgroup)
    tree.sort_descendants()
    return tree.write(format=9)


def processBPPresultPartA(filename, outgroup = 'J_neildi'):
    '''
    filename is a output file of BPP program
    return a list with trees and their posterior probability. Part A of the result.
    for the trees, in order to get consistent text representation of the trees, the outgroup is needed, and the nodes will be ordered
    
        (A) Best trees in the sample (3 distinct trees in all)
        74922 0.74922 0.74922 (D, (C, (A, B)));
        20016 0.20016 0.94938 (C, (D, (A, B)));
        5062 0.05062 1.00000 ((C, D), (A, B));
    with D as outgroup, return
        [ [0.74922, '(D,(C,(A,B)));'],
          [0.20016, '(D,((A,B),C));'],
          [0.05062, '(D,(C,(A,B)));']
        ]
    do not include trees with prob=0
    for trees with same topology, combine them
    '''
    ls_resultA = []
    f = open(filename)
    for line in f:
        if '(A) Best trees in the sample' in line:
            break
    
    for line in f:
        if line !='\n':
            ls_resultA.append(line)
        else:
            f.close()
            break
    
    resultA = []
    for line in ls_resultA:
        line = line.strip()
        counts, prob, cum_prob, tree = line.split(maxsplit=3)
        prob = float(prob)
        tree = orderTree(tree)
        if prob > 0:
             resultA.append([prob, tree])
    
    dc_result = {}
    for prob, tree in resultA:
        if tree not in dc_result:
            dc_result[tree] = []
        dc_result[tree].append(prob)
    for tree in dc_result:
        dc_result[tree] = sum(dc_result[tree])
    resultA = [[prob, tree] for tree, prob in dc_result.items()]
    resultA.sort(key=lambda x:x[0], reverse=True)
    
    return resultA


def combine2BPPresultPartA_List(resultA_bpp1, resultA_bpp2):
    '''
    resultA_bpp1, resultA_bpp2 are two list from processBPPresultPartA
    return a resultA by average the values
    '''
    dc_result = {}
    for prob, tree in resultA_bpp1 + resultA_bpp2:
        if tree not in dc_result:
            dc_result[tree] = []
        dc_result[tree].append(prob)
    for tree in dc_result:
        dc_result[tree] = sum(dc_result[tree]) / len(dc_result[tree])
    resultA = [[prob, tree] for tree, prob in dc_result.items()]
    resultA.sort(key=lambda x:x[0], reverse=True)
    return resultA

def combine2BPPresultPartAOne(bpp1, bpp2, outgroup = 'J_neildi', p_error = 0.3):
    '''
    bpp1 or bpp2 is files of bpp result, or two list of bpp result processed by processBPPresultPartA
    if the MAP species trees two files are the same, 
        select the two with highest MAP values. If their posterior probabilities differ <= p_error.
        if none-pairs have same MAP trees, if the mean absolute difference in the posterior probabilites is <= p_error, then combine them.
        if not, return None
    '''
    
    if type(bpp1) is str:
        resultA_bpp1 = processBPPresultPartA(bpp1, outgroup=outgroup)
    elif type(bpp1) is list:
        resultA_bpp1 = bpp1
    else:
        print('something wrong with',bpp1)
    if type(bpp2) is str:
        resultA_bpp2 = processBPPresultPartA(bpp2, outgroup=outgroup)
    elif type(bpp2) is list:
        resultA_bpp2 = bpp2
    else:
        print('something wrong with',bpp2)
    
    if len(resultA_bpp1) == 0 and len(resultA_bpp2) == 0:
        print(bpp1, bpp2, 'empty')
        return None
    if len(resultA_bpp1) == 0:
        print(bpp1,'empty')
        return resultA_bpp2
    if len(resultA_bpp2) == 0:
        print(bpp2,'empty')
        return resultA_bpp1
    
    if resultA_bpp1[0][1] == resultA_bpp2[0][1]:
        if abs(resultA_bpp1[0][0] - resultA_bpp2[0][0]) <= p_error:
            resultA = combine2BPPresultPartA_List(resultA_bpp1,resultA_bpp2)
            return resultA
    if resultA_bpp1[0][1] != resultA_bpp2[0][1]:
            tree1 = resultA_bpp1[0][1]
            tree2 = resultA_bpp2[0][1]
            dc_resultA1 = {tree:prob for prob, tree in resultA_bpp1}
            dc_resultA2 = {tree:prob for prob, tree in resultA_bpp2}
            if tree1 not in dc_resultA2:
                error1 = dc_resultA1[tree1]
            else:
                error1 = abs(dc_resultA1[tree1] - dc_resultA2[tree1])
            if tree2 not in dc_resultA1:
                error2 = dc_resultA2[tree2]
            else:
                error2 = abs(dc_resultA1[tree2] - dc_resultA2[tree2])
            if (error1 + error2)/2 < p_error:
                resultA = combine2BPPresultPartA_List(resultA_bpp1,resultA_bpp2)
                return resultA
    return None
        

def combineBPPresultPartAOne(files, outgroup = 'J_neildi', p_error = 0.3, outfile=None):
    '''
    files is a list with multiple filename of BPP result
    return a list of PartA by combining the results in these files
    there is at least two BPP results. If not, return None and print warning
    If there are more than two BPP results, choose two with highest maximum posterior probability (MAP) values.
        if the MAP species trees two files are the same, 
        select the two with highest MAP values. If their posterior probabilities differ <= p_error.
        if none-pairs have same MAP trees, if the mean absolute difference in the posterior probabilites is <= p_error, then combine them.
        else, use the other files until a good one is found.
        if not, return None and print warning
    if outfile is not None, write result to outfile, in the format of 'tree\tProb\n' in each line
    '''
    resultA_list = [processBPPresultPartA(f, outgroup=outgroup) for f in files]
    resultA_list = [e for e in resultA_list if len(e) > 0] #remove empty results
    fileNum = len(resultA_list)
    if fileNum == 0:
        print(files,'nothing in results')
        return None
    if fileNum == 1:
        return resultA_list[0]
    resultA_list.sort(key=lambda x:x[0][0], reverse=True) # sort results based on the best MAP values
    
    
    for j in range(1,fileNum):#try the pairs in the order of (0,1), (0,2),(1,2),(0,3),(1,3),(2,3)
        for i in range(j):
            resultA = combine2BPPresultPartAOne(resultA_list[i], resultA_list[j], outgroup = outgroup, p_error = p_error)
            if resultA is not None:
                if outfile is not None:
                    fout = open(outfile,'w')
                    for prob, tree in resultA:
                        fout.write('{tree}\t{prob:.6}\n'.format(tree=tree,prob=prob))
                    fout.close()
                return resultA
            
    print(files,'checked all but no good results')
    return None

