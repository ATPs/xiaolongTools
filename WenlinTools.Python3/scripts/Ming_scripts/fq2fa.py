import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

def find_dnlabel(alist):
    i = 0
    alist = [cmn.lastName(each) for each in alist]
    while (True):
        chars = [each[i] for each in alist]
        if len(set(chars)) == 1:
            i += 1
        else:
            break
    return alist[0][:i]
                

fq_list = sys.argv[1:]

dnlabel = find_dnlabel(fq_list).strip('_')
print(dnlabel)

dn = dnlabel + '_tmp'
dp = open(dn, 'w')

countfn = 0
for fq in fq_list:
    with open(fq) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 0:
                name = line.strip()
            elif i % 4 == 1:
                seq = line.strip()
                fasta = '>%s\n%s\n' % (name, seq)
                dp.write(fasta)
    countfn += i
dp.close()

#the maximal fasta of each fa
#if larger than this number, split them
#cutoff = 50000000
cutoff = 30000000 # entry

Nentry = (countfn + 1) / 4

Npack = ( Nentry / cutoff) + 1
    
#if Npack == 1:
if True:
    dnNew = dnlabel + '.fa'
    cmn.run('mv %s %s' % (dn, dnNew))
else:#need to split
    
    dnlist = []
    print('split files into %s packs' % Npack)
    each_pack = (Nentry / Npack + 1) * 2 #number of lines
    
    count = 1
    
    dnNew = '%s_%s.fa' % (dnlabel, count)
    dnlist.append(dnNew)
    dp = open(dnNew, 'w')
    with open(dn) as fp:
        for i, line in enumerate(fp):
            line = line.strip()
            if i % 2 == 0:
                defline = line
                continue
            elif i % 2 == 1:
                fasta = '%s\n%s\n' % (defline, line)
                dp.write(fasta)
            
            if (i+1) % each_pack == 0:
                count += 1
                dp.close()
                dnNew = '%s_%s.fa' % (dnlabel, count)
                dnlist.append(dnNew)
                dp = open(dnNew, 'w')
    
    dp.close()
    
    for each in dnlist:
        if not cmn.filexist(each):
            cmn.run('rm %s' % each)

    cmn.run('rm %s' % dn)
