#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_customBaits(customBaits):
    #for line in cmn.file2lines(fn):
    adict = {}
    for line in customBaits:
        sp, defline, seq = line.strip().split()
        #ensure primer
        if len(seq) != 658 and len(seq) != 698:
            print('Error! the length of the bait sequence is wrong for %s %s' % (sp, defline))
            sys.exit()

        if sp not in adict:
            adict[sp] = []
        if len(seq) == 698:
            adict[sp].append(line)
        else:# length is 658
            seq = add_primer(seq)
            adict[sp].append('%s\t%s\t%s' % (sp, defline, seq))
    return adict


def add_primer(seq):
    return 'ACTAATCATAAAGATATTGG%sTGATTTTTTGGTCATCCAGA' % seq.strip()


if __name__=='__main__':
    #options=parse_options()
    try:
        fqlist = sys.argv[1]
        fqlist = [os.path.abspath(each) for each in cmn.file2lines(fqlist)]
    except:
        print("Usage: *.py fqlist", file=sys.stderr)
        sys.exit()

    genusDict = {}
    try:
        for line in cmn.file2lines(sys.argv[2]):
            sp, genus = line.strip().split()
            if sp not in genusDict:
                genusDict[sp] = genus
    except:
        pass
    print(genusDict)
    #provided_baits = parse_customBaits(customBaits)

    fqDict = {}
    all_jobs = []
    for fq in fqlist:
        Id = cmn.lastName(fq).split('_')[0]
        Id = Id.replace('NVG-','').replace('11-BOA-','')
        try:
            fqDict[Id].append(fq)
        except:
            fqDict[Id] = [fq]

    cmn.run('rm */restricted_genus.info 2> /dev/null')
    finished = set([each.split('/')[-2].split('Run_')[-1] for each in cmn.cmd2lines('ls -d sampleRun_*/*ratio* 2> /dev/null')])
    for sample in fqDict:
        wdir = 'sampleRun_%s' % sample
        cmn.mkdir(wdir)
        os.chdir(wdir)

        #sampleInfo = sampleDict[sample] + '\n'
        #cmn.write_file(sampleInfo, 'sampleInfo')

        fqlist = fqDict[sample]
        cmn.write_lines(fqlist, 'fqlist')

        #try:
            #if baits are specified by the input
        #    baits = provided_baits[sample]
        #    cmn.write_lines(baits, 'sampleInfo.baits')
        #    cmd = 'bash /project/biophysics/Nick_lab/wli/sequencing/scripts/barcode_scripts/master_customBaits.sh %s' % sample
        #except KeyError:
        try:
            genus = genusDict[sample]
            cmn.write_file(genus, 'restricted_genus.info')
        except KeyError:
            pass

        if sample in finished:
            print('skip finished %s' % sample)
            continue

        cmd = 'bash /project/biophysics/Nick_lab/wli/sequencing/scripts/barcode_scripts/master.sh %s' % sample
        cmn.run('make_job.py "%s" -t 20 -p 256GB > bc%s.job' % (cmd, sample))
        #cmn.run('sbatch bc%s.job' % sample)
        os.chdir('..')
        all_jobs.append('%s/bc%s.job' % (wdir, sample))

    subCmds = ['sbatch %s\n' % job for job in all_jobs]
    cmn.write_file(''.join(subCmds), 'submission.sh')

    all_jobs.append('')
    cmn.write_lines(all_jobs, 'all_jobs')

