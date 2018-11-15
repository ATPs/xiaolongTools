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
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py IDlist", file=sys.stderr)
        sys.exit()


    exclusion = set(['All'])

    IDlist = set([each.replace('flt', '') for each in cmn.file2lines(fn)])

    takenList = set([])

    xl = pd.ExcelFile('/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/All-sequenced.xlsx')
    reference_table = '/project/biophysics/Nick_lab/wli/sequencing/scripts/downloaded_excels/reference_for_eachProject.txt'
    mapdict = {}
    for line in cmn.file2lines(reference_table):
        project, refs = line.strip().split()
        project = project.replace('_', ' ')
        mapdict[project] = refs.split(',')

    names = xl.sheet_names

    check_names = set(names) - exclusion

    sampleDict = {}
    for name in check_names:
        df = xl.parse(name)
        sheetIDs = list(df[df.columns[0]])
        sheetIDs = [str(each).replace('NVG-','').replace('11-BOA-', '').replace('LEP-', 'LEP')
                for each in sheetIDs]
        overlapping = set(sheetIDs) & IDlist
        if len(overlapping) != 0:
            for sample in overlapping:
                try:
                    sampleDict[sample].add(name)
                except KeyError:
                    sampleDict[sample] = set([name])
            takenList = takenList | overlapping

    #if 'Junonia' in todoList:
    #    print 'Note: please build Junonia coding regions.'

    #print 'the following tabs need to build trees:'
    #print '\n'.join(todoList)
    for sample in sampleDict:
        projects = sampleDict[sample]
        if 'Nymphalidae' in projects:
            projects.remove('Nymphalidae')
            if len(projects) == 0:
                print('Error: %s is only in Nymphalidae and would not be mapped' % sample, file=sys.stderr)
        for project in projects:
            #refs = mapdict[project]
            #for ref in refs:
            #    print sample, ref
            print(sample, project)
    missing = IDlist - takenList
    if len(missing) != 0:
        print('\n\nError! the following IDs can not be found in the table')
        print('\n'.join(missing))
