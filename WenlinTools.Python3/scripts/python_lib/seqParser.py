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
import numpy as np

class ST2(object):
    """
    This class parse ST2 result to a matrix
    """
    def __init__(self, fn):
        self.fn = fn
        self._parse_ST2()


    def _parse_ST2(self):
        lines = cmn.file2lines(self.fn)
        indices_start = False
        matrix_start = False
        self.pop2index = {}
        self.matrix = {}
        for line in lines:
            if 'Indices for populations:' in line:
                indices_start = True
                continue
            elif '----------------------' == line:
                indices_start = False

            if indices_start and ('---' not in line):
                #print line
                index, pop = line.strip().split()
                self.pop2index[pop] = index


            if 'Estimates for all loci ' in line:
                matrix_start = True
                continue
            elif matrix_start and line.strip() == '':
                matrix_start = False

            if matrix_start and ('===' not in line):
                items = line.strip().split()
                if items[0] == 'pop':
                    mheaders = items[1:]
                else:
                    i1 = items[0]
                    for ii, Fst in enumerate(items[1:]):
                        i2 = mheaders[ii]

                        key = self._make_key(i1, i2)
                        self.matrix[key] = float(Fst)

    def _make_key(self, i, j):
        alist = list(map(int, [i,j]))
        alist.sort()
        key = '-'.join(map(str, alist))
        return key

    def get_Fst_by_name(self, input_key):
        i, j = [self.pop2index[x] for x in input_key]

        Fst = self.get_Fst_by_index((i,j))
        return Fst

    def get_Fst_by_index(self, input_key):
        i, j = input_key

        key = self._make_key(i, j)
        Fst = self.matrix[key]
        return Fst


    def average_Fst(self, name):
        othernames = list(self.pop2index.keys())
        othernames.remove(name)

        Fst_list = [self.get_Fst_by_name((name, name2)) for name2 in othernames]
        return np.mean(Fst_list)

