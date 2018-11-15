import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

def parse_fqlist(fqlist):
    alist = []
    for line in cmn.file2lines(fqlist):
        label = cmn.lastName(line)
        if 'R1' in label:
            alist.append(line)
        elif 'R2' in label:
            alist.append(line)

    if len(alist) != 2:
        print('Error! can not recoginze fastq names in %s' % fqlist)
        cmd = 'touch fastq_error'
        cmn.run(cmd)
        sys.exit()
    return [alist]


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, outlabel, fqlist = sys.argv[1:]
    except:
        print("Usage: *.py reads.fa outlabel fqlist", file=sys.stderr)
        sys.exit()

    Ncpu = cmn.cpu_check()

    cmds = []
    cmd = 'platanus assemble -o %s -f %s -t %s -m 30' % (outlabel, fn, Ncpu)
    cmds.append(cmd)


    cmd = 'platanus scaffold -o %s -c %s_contig.fa -b %s_contigBubble.fa -t %s ' % (outlabel, outlabel, outlabel, Ncpu)

    #pairs = [each.split('---') for each in paired_libs.strip().split(',')]
    pairs = parse_fqlist(fqlist)

    for i, pair in enumerate(pairs):
        pair.sort()
        a, b = pair
        i = i + 1

        cmd += '-IP%s %s %s ' % (i, a, b)

    cmds.append(cmd)

    #do blast
    cmd = '\n\nmodule add blast'
    cmds.append(cmd)

    cmd = 'blastn -query %s_contig.fa -subject baits/bait0.fa -out bait0_denovo.br ' % (outlabel)
    cmd += ' -outfmt \'6 qseqid sseqid evalue pident length qlen qstart qend slen sstart send qseq sseq \''
    cmds.append(cmd)

    #check blast result
    cmd = 'python /project/biophysics/Nick_lab/wli/sequencing/scripts/barcode_scripts/check_blast_result.py bait0_denovo.br'
    cmds.append(cmd)

    cmds.append('\n')
    dn = 'denovo.cmds'
    cmn.write_lines(cmds, dn)

    cmd = 'bash denovo.cmds'
    cmn.run(cmd)
