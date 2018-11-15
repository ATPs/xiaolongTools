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
import os


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_pf_job(fns):
    cmds = []
    for fn in fns:
        cmd = "python /project/biophysics/Nick_lab/wli/sequencing/scripts/purity_filter.py %s > %s & " % (fn, cmn.lastName(fn))
        cmds.append(cmd)

    cmds.append('\nwait;\n')
    dn = 'pf.cmds'
    cmn.write_lines(cmds, dn)
    return dn

def make_mira_job(fns):
    cmds = []
    cwd = os.getcwd()
    new_dns = []
    for fn in fns:
        label = cmn.lastName(fn).replace('.fq', '')
        newdir = '%s/%s' % (cwd, label)
        cmd = 'mkdir %s;' % newdir
        cmd += 'cd %s;' % newdir
        cmd += '/home2/wli/local/mira_4.0.2_linux-gnu_x86_64_static/bin/mirabait -ik 30 -f fastq ../junk.fa ../%s %s >& %s.log &' % (fn, label, label)
        new_dns.append('%s/%s.fastq' % (newdir, label))
        cmds.append(cmd)
    cmds.append('\nwait;\n')

    os.chdir(cwd)
    dn = 'mira.cmds'
    cmn.write_lines(cmds, dn)
    return dn, new_dns


def make_qt_job(fns):
    cmds = []
    cwd = os.getcwd()
    new_dns = []
    for fn in fns:
        label = cmn.lastName(fn).replace('.fastq', '')
        cmd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/trim.py < %s > %s.fq &' % (fn, label)
        cmds.append(cmd)
        new_dns.append('%s/%s.fq' % (cwd, label))
    cmds.append('\nwait;\n')
    dn = 'qt.cmds'
    cmn.write_lines(cmds, dn)
    return dn, new_dns


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, label =sys.argv[1:3]
    except:
        print("Usage: *.py file_list_of_libs outdir_label", file=sys.stderr)
        sys.exit()

    ###############################################
    f_junk = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/junk.fa' #junk sequence filtered by mirabait
    ###############################################
    all_filenames = []

    #1. prepare data and directories
    wdir = 'preprocess_%s' % label
    cmn.mkdir(wdir)

    # link files
    fileNames = cmn.getid(fn)
    fns = [os.path.abspath(fileName) for fileName in fileNames]
    os.chdir(wdir)

    all_filenames.append(fns)

    lib_dir = '0_libs'
    cmn.mkdir(lib_dir)
    os.chdir(lib_dir)

    new_fns = []
    for fn in fns:
        cmn.run('ln -s %s' % fn)
        new_fns.append('../%s/%s' % (lib_dir, cmn.lastName(fn)))
    os.chdir('..')

    all_filenames.append(new_fns)
    #assuming demux has been done

    #1.purity filter
    pf_dir = '1_purify_filter'
    cmn.mkdir(pf_dir)
    os.chdir(pf_dir)

    job_file = make_pf_job(new_fns)
    cmn.run('bash %s' % job_file)
    os.chdir('..')

    new_dns = [i.replace(lib_dir, pf_dir) for i in new_fns]
    new_fns = new_dns
    all_filenames.append(new_fns)

    #2. mirabait
    mira_dir = '2_mirabait'
    cmn.mkdir(mira_dir)
    os.chdir(mira_dir)
    cmn.run('cp %s ./' % f_junk)

    #new_dns are parsed from the data
    job_file, new_dns = make_mira_job(new_fns)
    cmn.run('bash %s' % job_file)
    os.chdir('..')

    #new_dns = [i.replace(pf_dir, mira_dir) for i in new_fns]
    new_fns = new_dns
    all_filenames.append(new_fns)


    #3. quality trimmer
    qt_dir = '3_quality_trimmer'
    cmn.mkdir(qt_dir)
    os.chdir(qt_dir)

    #new_dns are parsed from the data
    job_file, new_dns = make_qt_job(new_fns)
    cmn.run('bash %s' % job_file)
    os.chdir('..')

    #new_dns = [i.replace(pf_dir, mira_dir) for i in new_fns]
    new_fns = new_dns
    all_filenames.append(new_fns)


    #4.re-pair IDs
    rp_dir = '4_re-pair'
    cmn.mkdir(rp_dir)
    os.chdir(rp_dir)

    #new_dns are parsed from the data
    job_file = 're-pair_auto.job'
    cmn.run('python /project/biophysics/Nick_lab/wli/sequencing/scripts/make_re-pair_cmds.py %s' % (' '.join(new_fns)))

    cmn.run('bash %s' % job_file)
    os.chdir('..')

    # fn change from label_R1.fq to label_paired1.fq
    new_dns = ['%s/%s' % (rp_dir, cmn.lastName(i).replace('R1', 'paired1')) for i in new_fns]
    new_fns = new_dns
    all_filenames.append(new_fns)

    cmn.pickle_write( all_filenames, 'filenames.list.pkl')
