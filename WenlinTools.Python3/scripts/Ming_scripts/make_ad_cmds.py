import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def make_ad_dict(fn):
    fi = '/work/biophysics/mtang/SNP_calling/scripts/data/adaptor/index_and_adaptor'
    #the index by the number
    adict = {}
    #the index by the barcode
    bdict = {}
    for line in cmn.file2lines(fi):
        index, barcode, seq = line.strip().split()
        adict[index] = seq
        bdict[barcode] = seq
    
    taken = set([])
    
    rdict = {}
    for line in cmn.file2lines(fn):
        try:
            sp, index = line.strip().split()
        except:
            continue

        if index.isdigit():
            ad1 = adict[index]
        else:
            try:
                ad1 = bdict[index]
            except KeyError:
                print('Error! can not find index info for %s' % line)
                sys.exit()

        rdict[sp] = ad1                

    return rdict
        
def group_fastqs(alist):
    adict = {}
    for each in alist:
        key = cmn.lastName(each).split('_')[0]
        try:
            adict[key].append(each)
        except KeyError:
            adict[key] = [each]
    
    #check and sort
    keys = list(adict.keys())
    for key in keys:
        each = adict[key]
        if len(each) != 2:
            print('Error! number of libs is wrong for %s' % key)
            print('below are the detected libs:')
            print('\n'.join(each))
            print('Please fix!')
            sys.exit()
        each.sort()
        adict[key] = each
    return adict        

#---------------main--------------------
fastqs = cmn.getid(sys.argv[1])

findex = sys.argv[2]

ad_dict = make_ad_dict(findex)

fastq_dict = group_fastqs(fastqs)

cmds = []
for key in fastq_dict:
    R1, R2 = fastq_dict[key]
    try:
        ad1 = ad_dict[key]
    except KeyError:
        print('Error! missing data for %s' % key)
        continue

    cmd = '/home2/wli/local/AdapterRemoval-1.5.4/AdapterRemoval --file1 %s --file2 %s --basename %s ' % (R1, R2, key)
    cmd += '--trimns --trimqualities --minquality 20 --pcr1 %s --pcr2 AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' % ad1
    cmds.append(cmd)

    #make job 
    fcmd = 'AD.cmds'
    cmn.write_lines(cmds, fcmd)
