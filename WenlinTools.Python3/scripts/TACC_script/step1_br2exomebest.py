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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def find_bestE_reads(fn, label):
    best_dict = {}
    with open(fn) as fp:
        for line in fp:
            items = line.strip().split()
            read_ID = items[1]
            evalue = float(items[3])
            line = '%s\t%s' % (label, line)
            try:
                current = best_dict[read_ID]
                if evalue < current[0]:
                    best_dict[read_ID] = (evalue, line)
            except:
                best_dict[read_ID] = (evalue, line)

    best_records = [i[1] for i in list(best_dict.values())]
    return best_records


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fns=sys.argv[1:]
    except:
        print("Usage: *.py filelist", file=sys.stderr)
        sys.exit()

    #1. get the best place for a hit
    #2. output the results by exome
    #3. put an taxon number infront of each line

    #outdir = 'blast_results'
    #cmn.run('rm -r %s' % outdir)
    #cmn.mkdir(outdir)

    #fns = cmn.getid(fn)
    print(fns)

    exome_dict = {}
    sp_dict = {}
    for fn in fns:
        label = cmn.lastName(fn).replace('.br', '').split('_')[0]
        print('parsing sample %s' % label)
        #best hit is the highest escore for each reads
        #also, add the label in the first place
        best_hits = find_bestE_reads(fn, label)

        #then separate it by exomes
        for line in best_hits:
            exon = line.strip().split()[1]#??
            try:
                exome_dict[exon].append(line)
            except:
                exome_dict[exon] = [line]

        try:
            sp_dict[label] += best_hits
        except:
            sp_dict[label] = best_hits

    print('writting outputs... ')
    cmn.pickle_write(exome_dict, 'blastByExon.dict.pkl')
    cmn.pickle_write(sp_dict, 'blastBySp.dict.pkl')
    #output the exome blasts
    #for exon in exome_dict:
    #    lines = exome_dict[exon]
    #    dn = '%s/%s.br' % (outdir, exon)
    #    cmn.write_file(''.join(lines), dn)







