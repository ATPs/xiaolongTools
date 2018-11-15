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



if __name__=='__main__':
    #options=parse_options()
    try:
        fqlist, label = sys.argv[1:]
    except:
        print("Usage: *.py fqlist label", file=sys.stderr)
        sys.exit()


    #pdir = '/archive/biophysics/Nick_lab/wli/workspace/proteinBaitBarocodes'

    fqs = cmn.file2lines(fqlist)

    #1. make blast db
    print('prepare libs...')
    db = '%s.fa' % label
    cmn.run('rm %s 2>/dev/null' % db)
    cmd = ''
    for fq in fqs:
        cmd += 'fq2fa %s >> %s;' % (fq, db)

    cmd += 'makeblastdb -in=%s -dbtype=nucl' % (db)
    cmn.run(cmd)


    print('running blast...')
    #2. do blast search
    Ncores = cmn.cpu_check()
    fout = label + '.cox1.out'
    cmd = "tblastn -task tblastn -max_target_seqs 50000000 -seg no -query /archive/biophysics/Nick_lab/wli/workspace/proteinBaitBarocodes/cox1.fa "
    cmd += "-db %s -out %s -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send evalue length pident nident qseq sseq' -num_threads %s" % (db, fout, Ncores)
    cmn.run(cmd)

    print('attempting assembling...')
    fID = label + '.blastID'
    cmd = 'cut -f 2 %s > %s' % (fout, fID)
    cmn.run(cmd)

    cmd = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/take_reps.py %s %s' % (db, fID)
    cmn.run(cmd)
    dn = label + '_taken.fa'

    cmd = 'platanus assemble -f %s -o %s.cox1 -t %s &> ass.log' % (dn, label, Ncores)
    cmn.run(cmd)

    #fass = '%s.cox1_contig.fa' % label
    #cmd = 'makeblastdb -in=%s -dbtype=nucl; tblastn -query /archive/biophysics/Nick_lab/wli/workspace/proteinBaitBarocodes/cox1.fa -db %s ' % (fass, fass)
    #cmd += '-out cox1_to_ass.out -outfmt \'6 qseqid sseqid qlen qstart qend slen sstart send\''
    #cmn.run(cmd)
