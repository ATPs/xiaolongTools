#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/biophysics/mtang/SNP_calling/script'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os
from collections import Counter


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


def separate_by_pair(fastqs):
    pdict = {}
    mapdict = { }
    for fastq in fastqs:
        key = '.'.join(cmn.lastName(fastq).split('.')[:-1])
        mapdict[key] = fastq

    
    names = list(mapdict.keys())
    length = max([len(name) for name in names])
    for i in range(length):
        if len(names) == 0:
            break

        checks = [name[:-1-i] for name in names]

        count_dict = Counter(checks)
        for key in count_dict:
            if count_dict[key] == 2: # got paired
                paired_names = [each for each in names 
                    if each.startswith(key)]
                fns = [mapdict[name] for name in paired_names]
                pdict[key] = fns

                #remove it from the list
                for name in paired_names:
                    names.remove(name)
            
    if len(pdict) == 0:
        print('Error! fastq lib name not recognized, contact Wenlin for help!')
        sys.exit()

    singleLibs = [mapdict[name] for name in names]
    print(singleLibs)
    
    if len(singleLibs) > 1:
        print('Warnning: more than one lib detected as single lib. below is the single list:')
        print('\n'.join(singleLibs))
        print('Email Wenlin for help')

    #print 'paired libs are:'
    #for key in pdict:
    #    print pdict[key]
    
    #print '\nsingle libs are:'
    #print ' '.join(singleLibs)

    singleFn = 'single.fq' 
    if cmn.filexist(singleFn):
        cmn.run('rm %s' % singleFn)
    for fn in singleLibs:
        cmn.run('cat %s >> %s' % (fn, singleFn))

    return pdict, singleFn



def separate_by_pair_old(label, fns):
    #paired = [i for i in fns if ('_paired' in i) ]
    paired = [i for i in fns if ('_pair' in i) or ('_R' in i)]
    paired.sort()
    if len(paired) != 2:
        print('error: wrong number of pairs as %s' % str(paired))
        print('from: %s' % str(fns))
        print('need to change the label criterion')
        sys.exit()

    unpaired = set(fns) - set(paired)

    #parse each files
    newPaired = []
    for fn in paired:
        if os.path.exists(cmn.lastName(fn)):
            cmn.run('unlink %s;' % cmn.lastName(fn))
        cmn.run('ln -s %s' % fn)
        newPaired.append(cmn.lastName(fn))

    singleFn = '%s_single.fq' % label
    if cmn.filexist(singleFn):
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
        finfo = sys.argv[1]
    except:
        print("Usage: *.py ../step1_gather_data/mapping_info.txt", file=sys.stderr)
        #print >> sys.stderr, "you should index assembly_v0.fa first with -p assembly_v0"
        #print >> sys.stderr, "using command /home2/wli/local/bwa-0.7.12/bwa index "
        sys.exit()

    cwd = os.getcwd()

    #1. read in info 
    refs = set([])
    rdict = {}
    for line in cmn.file2lines(finfo):
        sp, fastq, ref = line.strip().split()
        try:
            rdict[sp].append((fastq, ref))
        except KeyError:
            rdict[sp] = [(fastq, ref)]
        refs.add(ref)

    #2. prepare reference jobs
    refdir = '/work/biophysics/mtang/SNP_calling/indexed_references'
    cmn.mkdir(refdir)
    os.chdir(refdir)
    index_cmds = ['cd %s' % refdir]
    for ref in refs:
        if not os.path.exists(cmn.lastName(ref)):
            #cmn.run('ln -s %s' % ref)
            cmn.run('cp %s %s/' % (ref, refdir))
        ref = cmn.lastName(ref)
        reflabel = ref.replace('.fa', '')
        checkFn = reflabel + '.pac'
        if cmn.filexist(checkFn):
            print('found finished ref for %s, skip it' % ref)
            continue
        cmd = '/home2/wli/local/bwa-0.7.12/bwa index %s -p %s &' % (ref, reflabel)
        index_cmds.append(cmd)
    
    index_cmds.append('\nwait\n')
    os.chdir(cwd)
    
    print('#################################################')
    if len(index_cmds) != 2:
        dn = 'index.cmds'
        cmn.write_lines(index_cmds, dn)
        fjob = 'index.job'
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s -p 256GB > %s' % (dn, fjob)
        cmn.run(cmd)
        print('please submit %s to index references before start bwa mapping' % fjob)
    else:
        print('good news, all refs have been indexed')
    print('#################################################')

    #3. prepare mapping jobs

    cmn.mkdir('job_files')
    cmn.mkdir('cmd_files')

    todo_jobs = []
    for sp in rdict:
        records = rdict[sp]
        #print 'processing lib %s' % sp

        for record in records:
            fastq, ref = record
            reflabel = cmn.lastName(ref).replace('.fa', '')
            outlabel = '%s_%s' % (sp, reflabel)
            outdir = '%s/%s/%s' % (cwd, sp, reflabel)

            tmpcheck = cmn.cmd2lines(('ls %s/*sam 2> /dev/null' % outdir))
            if len(tmpcheck) > 0:
                total = 0
                for fn in tmpcheck:
                    total += cmn.filesize(fn)
                
                if total != 0:
                    print('skip finished mapping %s' % outdir)
                    continue

            cmn.mkdir(outdir)
            os.chdir(outdir)
        
            #paired is a dict
            paired, unpaired = separate_by_pair(fastq.split(','))

            cmd = 'cd %s;\n' % (refdir)
            for key in paired:
                lib1, lib2 = paired[key]
                cmd += '/home2/wli/local/bwa-0.7.12/bwa mem -t 32 -M %s %s %s > %s/%s_paired.sam;\n' % (reflabel, lib1, lib2, outdir, key)
            cmd += '/home2/wli/local/bwa-0.7.12/bwa mem -t 32 -M %s %s/%s > %s/single.sam;\n' % (reflabel, outdir, unpaired, outdir)
            #TODO: merge sams
            os.chdir(cwd)

            #tell which one is best
            #TODO: a script to tell who is best

            dn = 'cmd_files/%s.cmds' % outlabel
            cmn.write_file(cmd, dn)
            fjob = 'job_files/m%s.job' % outlabel
            cmn.run('/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s -t 20 > %s' % (dn, fjob))
            todo_jobs.append(fjob)
            
    
    print('#################################################')
    print('after finish index, please go to job_files directory and use the following commands to submit all jobs')
    print('')
    print('cd job_files')
    print('\n'.join(['sbatch %s' % cmn.lastName(fjob) for fjob in todo_jobs]))
        #cmn.run('/home2/wli/my_programs/decorate_job.py %s -p 256GB > %s' % (dn, fjob))
        #cmn.run('cd job_files; sbatch %s.job' % label)






