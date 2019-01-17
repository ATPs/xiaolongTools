

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

def extractStructureTableFromCommandlineOutput(filename, output=None, processLineForR=True):
    '''
    filename is the output run from STRUCTURE commandline
    
Inferred ancestry of individuals:
        Label (%Miss) :  Inferred clusters
  1     5490    (3)   :  0.998 0.002 
  2 16106C04    (9)   :  0.906 0.094 
  3     3935    (0)   :  1.000 0.000 
  4     4758    (1)   :  1.000 0.000 
    if processLineForR is True
    extract the table as
    Label    K2_1    K2_2
    5490     0.906   0.002
    
    also return a list of info [Label, K2_1, K2_2]
    else return the line directly
    1     5490    (3)   :  0.998 0.002 
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
    if processLineForR:
        results = list(map(processSTRUCTUREResultLine, ls_info))
    else:
        results = ls_info
    
    if output is not None:
        fout = open(output,'w')
        for e in results:
            if processLineForR:
                fout.write('\t'.join(str(i) for i in e)+'\n')
            else:
                fout.write(e)
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
    parser.add_argument('-r','--processLineForR',help = 'whether to change the lines so that it can be used by R or python directly. default=1. (0 for False, 1 for True)', default = 1, choices=[0,1], type=bool)
    f = parser.parse_args()
    extractStructureTableFromCommandlineOutput(filename=f.input, output=f.output,processLineForR=f.processLineForR)
