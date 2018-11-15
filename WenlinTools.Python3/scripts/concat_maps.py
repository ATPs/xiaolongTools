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


def read_length_info(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        scaf, length = line.strip().split()
        adict[scaf] = int(length)
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_index_header(fvcf, flen):

    lendict = read_length_info(flen)

    ordered_scafs = []
    taken = set([])
    with open(fvcf) as fp:
        for line in fp:
            if not line.startswith('scaffold'):
                continue
            scaf = line.split()[0]
            if scaf not in taken:
                taken.add(scaf)
                ordered_scafs.append(scaf)

    info = []
    for scaf in ordered_scafs:
        length = lendict[scaf]
        for i in range(length):
            info.append('%s\t%s\n' % (scaf, (i+1)))

    dn = 'index_header'
    cmn.write_file(''.join(info), dn)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    fns = cmn.cmd2lines('ls *.map| grep -v all| grep -v test| grep -v concat')
    #fns = fns[:1]

    f_label = '15101E05_snp.vcf'


    #make the scaffold index
    #cmd = "grep -v '^#' 15101E05_snp.vcf| cut -f 1,2 > index_header"
    if not cmn.filexist('index_header'):
    #    cmn.run(cmd)
        make_index_header('15101E05_snp.vcf', 'assembly_v2_length.txt')

    header = ['scaffold', 'index']

    for fn in fns:
        label = fn.split('_')[0]
        header.append(label+'_f')
        header.append(label+'_m')

    cmn.write_lines(fns, 'map_name_order')

    cmn.write_file('\t'.join(header)+'\n', 'table_header')

    cmd = 'cp table_header all_concat.map;'
    cmd += 'paste index_header %s >> all_concat.map;' % ' '.join(fns)
    if not cmn.filexist('all_concat.map'):
        print(cmd)
        cmn.run(cmd)
    else:
        print('the final file all_concat.map has exist!, skip!')
