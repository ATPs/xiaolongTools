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


def separate_by_pair(fastqs, wdir):
    print(wdir)
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

    singleFn = '%s/single.fq' % wdir
    if cmn.filexist(singleFn):
        cmn.run('rm %s' % singleFn)
    for fn in singleLibs:
        cmn.run('cat %s/%s >> %s' % (wdir, fn, singleFn))

    return pdict, singleFn


def separate_by_pair_vold(fastqs, wdir):
    pdict = {}
    mapdict = { }
    for fastq in fastqs:
        key = '.'.join(cmn.lastName(fastq).split('.')[:-1])
        mapdict[key] = fastq

    
    names = list(mapdict.keys())
    length = min([len(name) for name in names])
    for i in range(length):
        checks = [name[:-1-i] for name in names]
        count_dict = Counter(checks)
        if max(count_dict.values()) == 2: #got paired
            for name in names:
                key = name[:-1-i]
                fn = mapdict[name]
                try:
                    pdict[key].append(fn)
                except:
                    pdict[key] = [fn]
            break
    
    if len(pdict) == 0:
        print('Error! fastq lib name not recognized, contact Wenlin for help!')
        sys.exit()

    singleLibs = []
    keys = list(pdict.keys())
    for key in keys:
        libs = pdict[key]
        if len(libs) != 2:
            singleLibs += libs
            del pdict[key]
                
    #print 'paired libs are:'
    #for key in pdict:
    #    print pdict[key]
    
    #print '\nsingle libs are:'
    #print ' '.join(singleLibs)

    singleFn = '%s/single.fq' % wdir
    if cmn.filexist(singleFn):
        cmn.run('rm %s' % singleFn)
    for fn in singleLibs:
        cmn.run('cat %s >> %s' % (fn, singleFn))

    return pdict, singleFn


def make_bwa_cmds(fastqs, reflabel, wdir):
    outdir = '.'
    #1. index it
    #this has been done in the template

    #2. pairing dict
    paired, unpaired = separate_by_pair(fastqs, wdir)

    fsams = []
    cmd = ''
    for key in paired:
        lib1, lib2 = paired[key]
        fsam = '%s/%s_paired.sam' % (outdir, key)
        cmd += 'bwa mem -t 16 -M %s %s %s > %s;\n' % (reflabel, lib1, lib2, fsam)
        fsams.append(fsam)

    cmd += 'bwa mem -t 16 -M %s %s/%s > %s/single.sam;\n' % (reflabel, outdir, cmn.lastName(unpaired), outdir)
    fsams.append('%s/single.sam' % outdir)

    return cmd, fsams


def tell_best_ref(fmap, refs):
    aset = set(refs)
    final = (None, 0)
    for line in cmn.file2lines(fmap):
        items = line.strip().split()
        ref = items[0]
        mapped = int(items[1])

        if ref not in aset:
            continue

        if mapped > final[1]:
            final = (ref, mapped)


    if final[1] == 0:
        print('Error! no info to tell which is best for %s' % str(aset))
        sys.exit()

    return final[0]    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        freftable, mapdir, freq, fsubset = sys.argv[1:5]
    except:
        print("Usage: *.py ../step1_gather_data/mapping_info.txt ../step2_bwa_mapping ../step1_gather_data/require_SNPs.dict.pkl TACC_IDs", file=sys.stderr)
        sys.exit()

    cwd = os.getcwd()

    #subsetIDs = set(cmn.getid(fsubset))
    subsetJobs = set([cmn.lastName(line.replace('sbatch', '').strip())[4:-4]
                for line in cmn.file2lines(fsubset)])

    #1. read in info
    fsams = cmn.cmd2lines('ls %s/*/*/*.sam' % mapdir)
    #print fsams
    samdirs = set(['/'.join(fsam.split('/')[:-2]) for fsam in fsams])
    #print samdirs
    require_refs = cmn.pickle_read(freq)

    fq_dict = {}
    refdict = {}
    #1. tell by reftable
    #make the requirement by the reftable
    for line in cmn.file2lines(freftable):
        items = line.strip().split()
        sp = items[0]
        fastqs = items[1].split(',')
        fq_dict[sp] = fastqs

    # check for reference
    #2. tell by best mapping
    for samdir in samdirs:
        sp = cmn.lastName(samdir)
        if sp not in refdict:
            refdict[sp] = []
        require_list = require_refs[sp]
        
        for refs in require_list:
            refs = ['.'.join(cmn.lastName(ref).split('.')[:-1]) for ref in refs]

            if len(refs) == 1:
                ref = refs[0]
                refdict[sp].append((samdir, ref))
            else:
                fmap = '%s/mapping_stat.txt' % samdir
                
                #fbest = '%s/best_mapping.txt' % samdir
                #if not os.path.exists(fbest):
                if not os.path.exists(fmap):
                    print('checking ref for %s' % sp)
                    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/tell_best_mapping.py %s' % samdir
                    #print cmd
                    cmn.run(cmd)
                ref = tell_best_ref(fmap, refs)

                refdict[sp].append((samdir, ref))
    # if the ref info contain the ref in multiple lines, use it
    
    #2. make gatk index of reference
    
    #this ref dict is universal to avoid duplicated indexes
    refdir = '/work/biophysics/mtang/SNP_calling/indexed_references'
    
    cmds = ['cd %s' % refdir]
    cmds.append('module add picard/1.117')
    cmds.append('module load java/oracle/jdk1.8.0_65')
    
    todoref_count = 0
    taken_refs = set([])
    for sublist in list(refdict.values()):
        for samdir, ref in sublist:
            if ref in taken_refs:
                continue
            else:
                taken_refs.add(ref)
            fcheck = '%s/%s.dict' % (refdir, ref)
            #print refdir, ref
            if cmn.filexist(fcheck):
                print('skip finished indexed %s' % ref)
                continue
            # this has been finished in bwa mapping    
            #/home2/wli/local/bwa-0.7.12/bwa index -p assembly_selfref assembly_selfref.fa > index.log &
            todoref_count += 1
            cmds.append('java -jar $PICARD/CreateSequenceDictionary.jar R=%s.fa O=%s.dict &' % (ref, ref))

            cmds.append('/home2/wli/local/samtools-1.2/samtools faidx %s.fa &' % ref)

    cmds.append('\nwait;\n')

    isIndexed = False
    print('###############################################')
    if todoref_count != 0:
        fcmd = 'gatkIndex.cmd'
        cmn.write_lines(cmds, fcmd)
        fjob = 'gatkIndex.job'
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s -p 256GB > %s' % (fcmd, fjob)
        cmn.run(cmd)
        print('please submit %s to the queue for indexing ' % fjob)
    else:
        print('good news! all references have been indexed') 
        isIndexed = True
    print('###############################################')
    
    if not isIndexed:
        print('**********************************************')
        print('\nimportant!!!')
        print('please re-run this script after all references are indexed!\n')
        print('**********************************************')
    ###############################
    #all the steps below would put into the job files

    template = cmn.txt_read('/work/biophysics/mtang/SNP_calling/scripts/templates/template_gatk_unbias4TACC.job')
    

    cmn.mkdir('job_files')
    fjobs = []
    for sp in refdict:
        #if sp.split('_')[0] not in subsetIDs:
        #    continue
        snp_list = refdict[sp]
        for samdir, ref in snp_list:
            label = '%s_%s' % (sp, ref)
            if label not in subsetJobs:
                continue
            print('processing %s' % label)
            
            #a. make directory
            olabel = '%s_%s' % (sp, ref)
            wdir = '%s/%s' % (cwd, olabel)
            wdir4TACC = '../%s' % olabel
            cmn.mkdir(wdir)

            #samdir, ref = refdict[sp]
            print(samdir, ref)
            fsams = cmn.cmd2lines('ls %s/%s/*.sam' % (samdir, ref))
            #for tacc: copy sam files into this directory
            print('copying sam files for %s...' % olabel)
            for fsam in fsams:
                cmd = 'cp %s %s' % (fsam, wdir)
                cmn.run(cmd)
            fsams = [cmn.lastName(fsam) for fsam in fsams]    

            #for tacc: copy all the referenece info into this directory
            print('copying reference info for %s...' % olabel)
            cmd = 'cp %s/%s* %s' % (refdir, ref, wdir)
            cmn.run(cmd)

            fass = '%s.fa' % (ref)
            #wdir = '%s/%s_%s' % (cwd, sp, ref)
        
            info = template.replace('5642', olabel)
            info = info.replace('[WL_ref]', fass)
        
    #3. preparing the directories
            cmds = ['cd %s\n\n' % wdir4TACC]

        #merge sam files into one
            cmd = 'python /work/00412/mtang/sequencing/scripts/merge_mapped_sams.py %s.sam %s' % (olabel, ' '.join(fsams))
            cmds.append(cmd)
        
            cmds.append('\ncd %s\n' % (wdir4TACC))

            info = info.replace('[WL_preprossing]', '\n'.join(cmds))
    
    #4. first run of gatk to the realign.bam
    #this step has been fully included in the template

    #5. correct bias
    #this step has been fully included in the template

    #6. realign sams to the new reference genome
    #the new reference is named 'assembly_selfref_v2.fa'
            fastqs = fq_dict[sp]
            #for tacc, copy fastq to wdir
            print('copying fastq for %s...' % olabel)
            for fastq in fastqs:
                cmd = 'cp %s %s' % (fastq, wdir)
                cmn.run(cmd)
            
            fastqs = [cmn.lastName(fastq) for fastq in fastqs]
            cmd, fnewSams = make_bwa_cmds(fastqs, 'assembly_selfref_v2', wdir)
            info = info.replace('[WL_mapping_cmds]', cmd)
    #7. merge mapped sams
            cmd = 'python /work/00412/mtang/sequencing/scripts/merge_mapped_sams.py %s_step2.sam %s' % (olabel, ' '.join(fnewSams))
            info = info.replace('[WL_merge_sam_cmds]', cmd)

    #8. re-run gatk
    #this step has been fully included in the template
            fjob = 'job_files/gatk%s.job' % olabel
            cmn.write_file(info, fjob)

            fjobs.append(fjob)


    print('##########################################################################')        
    print('please use the following cmd to submit unfinished jobs')
    print('cd job_files')
    print('\n'.join(['sbatch %s' % cmn.lastName(fjob) for fjob in fjobs]))
    print('##########################################################################')        

