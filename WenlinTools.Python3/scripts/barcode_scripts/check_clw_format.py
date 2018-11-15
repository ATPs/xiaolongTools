#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    lenDict = {}
    for line in cmn.file2lines(fn):
        if '(assembled)' in line:
            name, seq = line.strip().split()
            if '_' not in name:
                print('did\'t change name for %s' % name)

            if seq.count('N') != 0 or seq.count('-') != 0 or seq.count('X') != 0:
                print('missing letter in %s' % name)

            if not seq.isupper():
                print('small letter exists in %s' % name)
            length = len(seq)

            try:
                lenDict[length].append(name)
            except KeyError:
                lenDict[length] = [name]
        if '_SRNP_' in line:
            print('warning! change SRNP number for this line: %s' % line)
    if len(lenDict) != 1:
        for key in lenDict:
            print(key, lenDict[key])




