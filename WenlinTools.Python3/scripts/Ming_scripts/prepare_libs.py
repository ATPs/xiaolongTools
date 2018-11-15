import sys, os

sys.path.append('/work/biophysics/mtang/SNP_calling/scripts')

import cmn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_info_file(wdir):
    fns = cmn.cmd2lines(' ls %s/fastq*' % wdir)
    adict = {}
    for fn in fns:
        fp = open(fn)
        #with open(fn) as fp:
        for line in fp:
            line = line.strip()
            sp = cmn.lastName(line).split('_')[0]
            try:
                adict[sp].append(line)
            except KeyError:
                adict[sp] = [line]
        fp.close()       

    for key in list(adict.keys()):
        adict[key] = set(adict[key])

    return adict


def parse_fastq_dir(wdir):
    adict = {}
    fns = cmn.cmd2lines('ls %s/*q' % wdir)
    for line in fns:
        sp = cmn.lastName(line).split('_')[0]
        try:
            adict[sp].append(line)
        except KeyError:
            adict[sp] = [line]
    return adict



def parse_done_file(fn):
    adict = {}
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            sp = line.split()[0].split('_')[0]
            try:
                adict[sp].append(line)
            except KeyError:
                adict[sp] = [line]
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


reftable = sys.argv[1]
fastqs = [line.strip() for line in cmn.file2lines(sys.argv[2])]

#1. read in refdict
#1.2 also get the requirement for each sample
refdict = {}
allrefs = set([])
requires = {}
for line in cmn.file2lines(reftable):
    if line[0] == '#':
        continue
    if line.strip() == '':
        continue
    items = line.strip().split()
    sp = items[0]
    refs = items[1:]
    allrefs = allrefs | set(refs)
    try:
        refdict[sp] += refs
    except KeyError:
        refdict[sp] = refs

    try:
        requires[sp].append(set(refs))
    except KeyError:
        requires[sp] = [set(refs)]

#check if the ref genome exist
#check if the ref is conflict with the one we already have
refdir = '/work/biophysics/mtang/SNP_calling/indexed_references'
for ref in allrefs:
    if not os.path.exists(ref):
        print('reference %s doesn\'t exist! please email to ask!' % ref)
    
    oldref = '%s/%s' % (refdir, cmn.lastName(ref))
    #print oldref, ref
    if os.path.exists(oldref):
        check = cmn.cmd2info('diff %s %s| wc -l ' % (oldref, ref))
        if int(check) != 0:
            print('new ref is different from old ref! please email to ask!')
            print('old ref: %s' % oldref)
            print('new ref: %s' % ref)

#addon: check fastq to see if anything has been done before
fdone = '/project/biophysics/Nick_lab/mtang/archive/submission_done'
info_wdir = '/project/biophysics/Nick_lab/mtang/archive/step1_info'
fastq_dir = '/project/biophysics/Nick_lab/mtang/archive/fastq_libs'

info_dict = parse_info_file(info_wdir)
done_dict = parse_done_file(fdone)
#current_fastqs = parse_fastq_dir(fastq_dir)

print('***********************************************************')
print('please email the following info to whoever send you request\n\n')
sp_list = set([cmn.lastName(fastq).split('_')[0] for fastq in fastqs])

for sp in sp_list:
    #currently, only use sp to check
    #check the info file
    hasDone = False
    if sp in done_dict:
        print('sample %s has been done before in Qian\'s spreedsheet, here is the info:' % sp)
        print('\n'.join(done_dict[sp]))
        hasDone = True

    #check the previous run info 
    if sp in info_dict:
        print('sample %s has been run before, here is the info from previous run:' % sp)
        print('\n'.join(info_dict[sp]))
        hasDone = True

    #check current fastq
    #if sp in current_fastqs:
    #    print 'found fastq for sample %s, please make sure you combined data:' % sp
    #    print '\n'.join(current_fastqs[sp])
    #    hasDone = True

    if hasDone:
        print('Please make sure you have combined data for %s' % sp)
    else:
        print('Do not find previous data for %s, did we miss it?' % sp)
    print('\n')

print('***********************************************************')
    


#2. group fastq by SP
qdict = {}
for fastq in fastqs:
    if not os.path.exists(fastq):
        print('fastq file %s doesn\'t exist! please email to ask!' % fastq)
    sp = cmn.lastName(fastq).split('_')[0]
    try:
        qdict[sp].append(fastq)
    except KeyError:
        qdict[sp] = [fastq]

#3. check sp to see if refs are specified
goodSPs = []
for sp in qdict:
    if sp in refdict:
        goodSPs.append(sp)
    else:
        print('no reference genome found for sample %s, please email to ask' % sp, file=sys.stderr)    

#4. output the mapping relationship
new = []
for sp in goodSPs:
    for ref in refdict[sp]:
        new.append('%s\t%s\t%s\n' % (sp, ','.join(qdict[sp]), ref))

dn = 'mapping_info.txt'
cmn.write_file(''.join(new), dn)


dn = 'require_SNPs.dict.pkl'
cmn.pickle_write(requires, dn)
