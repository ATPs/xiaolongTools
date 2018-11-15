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
from fullname_lib import get_names


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def find_SRNPnumber(line):
    for item in line.strip().split('|'):
        if "SRNP" in item:
            return item


def update_SRNP_species(line2):
    srnp = find_SRNPnumber(line2)
    cmd = '/home2/wli/anaconda/bin/python /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/updateSRNPnumber.py %s' % srnp
    print(cmd)
    sp = cmn.cmd2info(cmd).strip()
    if sp != '':
        newname = '%s|%s' % ('_'.join(sp.split()), srnp)
    else:
        newname = line2
    return newname


if __name__=='__main__':
    #options=parse_options()
    try:
        sampleID = sys.argv[1]
    except:
        print("Usage: *.py sampleID", file=sys.stderr)
        sys.exit()

    cmd = 'verify_barcodes_4checking.py %s_exactCheck.fasta > verify.log' % sampleID
    cmn.run(cmd)


    fa = '%s_exactCheck.fasta_paired.fa' % sampleID
    freport = '%s_exactCheck.fasta.report' % sampleID

    lines = cmn.cmd2lines('grep ">" %s' % fa)
    lines1 = []
    lines2 = []
    rdict = {}
    for i, line in enumerate(lines):
        line = line[1:]
        if 'assembled' in line:
            lines1.append(line)
            justTake = True
        else:
            if justTake == True:
                lines2.append(line)
                justTake = False
            else:
                lines2[-1] = line

    diff_chars = {}
    for line in cmn.file2lines(freport):
        items = line.split()
        #sp = items[0].split('_')[0]
        sp = items[0]
        try:
            diff_chars[sp] += 1
        except KeyError:
            diff_chars[sp] = 1


    nameDict = get_names()
    #for line in cmn.file2lines('IDlist.renamed'):
    #    a, b = line.strip().split()
    #    nameDict[a.split('_')[0]] = b

    dn = '%s_summaryCheck.table' % sampleID
    dp = open(dn, 'w')
    dp.write('sample\tcheckFlag\tdiffN\tbarcodeSP\tannotatedSP\n')
    for line1, line2 in zip(lines1, lines2):
        #print 'line1', line1
        #print 'line2', line2
        sp = line1.replace('(assembled)', '')
        label = sp.split('_')[0]
        if line2[:5] == 'JTRIO':
            line2 = 'Calephelis_%s' % line2
        if 'SRNP' in line2:
            line2 = update_SRNP_species(line2)

        if '_' not in line2 and '|' not in line2:
            line2 = line2.replace('-', '_')
        items = line2.split('_')

        if any([char.isdigit() for char in items[0]]):
            ref = items[1:3]
        else:
            ref = items[:2]

        if len(ref) == 0:
            ref = items[0]
        else:
            genus = ref[0]
            ref = '_'.join(ref).split('|')[0]
        try:
            Ndiff = diff_chars[sp]
        except KeyError:
            Ndiff = 0

        if sp in rdict:
            cur = rdict[sp]
            if cur[2] > Ndiff:
                rdict[sp] = [ (genus in nameDict[label]), Ndiff, line2, nameDict[label].replace(' ', '_').replace('\t', '_') ]
        else:
            rdict[sp] = [ (genus in nameDict[label]), Ndiff, line2, nameDict[label].replace(' ', '_').replace('\t', '_') ]


    for sp in rdict:
        dp.write('%s\t%s\n' % (sp, '\t'.join(map(str, rdict[sp]))))
    dp.close()

