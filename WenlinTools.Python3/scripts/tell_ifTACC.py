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
def parse_refTable(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        items = line.strip().split()
        sample = items[0]
        refs = items[1:]
        maxSize = 0
        for fref in refs:
            size = compute_fileSize([fref])
            maxSize = max(size, maxSize)

        adict[sample] = maxSize
    return adict

def compute_fileSize(alist):
    size = 0
    for fn in alist:
        if 'archive/butterfly' in fn:
            cmd = 'ssh wenlin@alea.swmed.edu "python /home/wenlin/my_programs/filesize.py %s"' % fn
            size += int(cmn.cmd2info(cmd).strip())
        else:
            size += cmn.filesize(fn) / 1024 / 1024
    return size

def parse_fastqlist(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        sample = cmn.lastName(line).split('_')[0]
        try:
            adict[sample].append(line)
        except KeyError:
            adict[sample] = [line]
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fref, fqlist = sys.argv[1:3]
    except:
        print("Usage: *.py reftable.txt fastq.filelist", file=sys.stderr)
        sys.exit()


    refCut = 500 # MB
    fqCut = 15000 # 15GB total

    refSize_dict = parse_refTable(fref)

    fq_dict = parse_fastqlist(fqlist)

    for sample in fq_dict:
        fqs = fq_dict[sample]

        fqSize = compute_fileSize(fqs)
        refSize = refSize_dict[sample]

        if refSize <= refCut and fqSize <= fqCut:
            label = 'TACC'
        else:
            if refSize <= 650 and fqSize <= 5000:
                label = 'TACC'
            else:
                if fqSize <= 8000:
                    label = 'stampede2'
                else:
                    label = 'bioHPC'

        print(sample, '%.1fGB' % (fqSize/1024), '%.1fMB' % refSize, label)




