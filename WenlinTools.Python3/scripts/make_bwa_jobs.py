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
def separate_by_label(fns):
    adict = {}
    for fn in fns:
        label = '_'.join(cmn.lastName(fn).split('_')[:-1])
        try:
            adict[label].append(fn)
        except:
            adict[label] = [fn]
    return adict


def separate_by_pair(label, fns):
    paired = [i for i in fns if ('_paired' in i) or ('_R' in i)]
    paired.sort()
    if len(paired) != 2:
        print('error: wrong number of pairs as %s' % str(paired))
        print('from: %s' % str(fns))
        #print 'need to change the label criterion'
        print('skip this lib')
        return None, None

    unpaired = set(fns) - set(paired)

    #parse each files
    newPaired = []
    for fn in paired:
        cmn.run('ln -s %s' % fn)
        newPaired.append(cmn.lastName(fn))

    singleFn = '%s_single.fq' % label
    if cmn.filexist:
        cmn.run('rm %s' % singleFn)
    for fn in unpaired:
        cmn.run('cat %s >> %s' % (fn, singleFn))

    return newPaired, singleFn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        odir, f_ass = sys.argv[1:3]
    except:
        print("Usage: *.py filelist assembly_v0.fa", file=sys.stderr)
        print("you should index assembly_v0.fa first with -p assembly_v0", file=sys.stderr)
        print("using command /home2/wli/local/bwa-0.7.12/bwa index ", file=sys.stderr)
        sys.exit()

    #fns = cmn.cmd2lines('ls %s/*.fq' % odir)
    fns = cmn.getid(odir)

    group_dict = separate_by_label(fns)

    ass_label = cmn.find_between(cmn.lastName(f_ass), 'assembly_', '.fa')

    cmn.mkdir('job_files')
    cmn.mkdir('cmd_files')

    for plabel in group_dict:
        print('processing lib %s' % plabel)
        each = group_dict[plabel]
        #also parse the files inside this function
        #return the file name after parsing
        paired, unpaired = separate_by_pair(plabel, each)
        if paired == None:
            continue
        label = '%s_%s' % (plabel, ass_label)
        #index_label = cmn.lastName(f_ass).replace('.fa', '')
        index_label = f_ass.replace('.fa', '')
        cmd = ''
        cmd += '/home2/wli/local/bwa-0.7.12/bwa mem -t 32 -M %s %s %s > %s_paired.sam;\n' % (index_label, paired[0], paired[1], label)
        cmd += '/home2/wli/local/bwa-0.7.12/bwa mem -t 32 -M %s %s > %s_single.sam;\n' % (index_label, unpaired, label)
        dn = 'cmd_files/%s.cmds' % label
        cmn.write_file(cmd, dn)
        fjob = 'job_files/m%s.job' % label
        cmn.run('/home2/wli/my_programs/decorate_job.py %s > %s' % (dn, fjob))






