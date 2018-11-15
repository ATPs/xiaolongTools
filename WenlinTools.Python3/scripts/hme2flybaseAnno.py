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

def find_annotation(CG):
    if CG == 'NA':
        return 'NA'
    fn = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/gene_association.fb'
    anno_lines = cmn.cmd2lines('grep %s %s' % (CG, fn))
    annos = set([])
    for line in anno_lines:
        items = line.split('\t')
        anno = items[9]
        annos.add(anno)
    return ','.join(annos)


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py IDlisti/countlist", file=sys.stderr)
        sys.exit()

    lines = cmn.file2lines(fn)

    rawhmeIDs = [line.split()[-1] for line in lines]
    #expect ID as HMEL015065-RA
    hmeIDs = [each.split('-')[0][4:] for each in rawhmeIDs]

    hme2CG = {}
    fn = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/hme_to_flybase'
    for line in cmn.file2lines(fn):
        #hme000003	CG15845-PE
        hme, CG = line.strip().split()
        hme = hme.replace('hme', '')
        CG = CG.split('-')[0]
        hme2CG[hme] = CG

    CG_list = []
    for ID in hmeIDs:
        try:
            CG = hme2CG[ID]
        except:
            CG = 'NA'
        CG_list.append(CG)


    for i, CG in enumerate(CG_list):
        anno = find_annotation(CG)
        info = '%s\t%s\t%s' % (rawhmeIDs[i], CG, anno)
        newline = lines[i].replace(rawhmeIDs[i], info)
        print(newline)




