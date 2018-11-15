

def processSTRUCTUREResultLine(line):
    '''
    line = ' 67     4507    (1)   :  1.000 0.000 \n'
    return ['4507',1,0]
    '''
    l,r = line.split(':')
    name = l.split()[1]
    values = r.split()
    values = [float(e) for e in values]
    return [name] + values

def extractStructureTableFromCommandlineOutput(filename, output=None):
    '''
    filename is the output run from STRUCTURE commandline
    
Inferred ancestry of individuals:
        Label (%Miss) :  Inferred clusters
  1     5490    (3)   :  0.998 0.002 
  2 16106C04    (9)   :  0.906 0.094 
  3     3935    (0)   :  1.000 0.000 
  4     4758    (1)   :  1.000 0.000 

    extract the table as
    Label    K2_1    K2_2
    5490     0.906   0.002
    
    also return a list of info [Label, K2_1, K2_2]
    '''
    f = open(filename)
    for line in f:
        print(line)
        if 'Inferred ancestry of individuals:' in line:
            break
    f.readline()
    ls_info = []
    for line in f:
        if line != '\n':
            ls_info.append(line)
        else:
            f.close()
            break
    
    results = map(processSTRUCTUREResultLine, ls_info)
    
    if output is not None:
        fout = open(output,'w')
        for e in results:
            fout.write('\t'.join(str(i) for i in e)+'\n')
        fout.close()
    
    return results


description = '''filename is the output run from STRUCTURE commandline
    
Inferred ancestry of individuals:
        Label (%Miss) :  Inferred clusters
  1     5490    (3)   :  0.998 0.002 
  2 16106C04    (9)   :  0.906 0.094 
  3     3935    (0)   :  1.000 0.000 
  4     4758    (1)   :  1.000 0.000 

    extract the table as
    Label    K2_1    K2_2
    5490     0.906   0.002
    
    also return a list of info [Label, K2_1, K2_2]
'''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file of STRUCTURE output', required=True)
    parser.add_argument('-o','--output',help = 'where to store the output. default=None, only return the result', default = None)
    f = parser.parse_args()
    extractStructureTableFromCommandlineOutput(filename=f.input, output=f.output)
