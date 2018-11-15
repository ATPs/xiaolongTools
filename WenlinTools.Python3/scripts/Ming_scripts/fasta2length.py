import sys

fn = sys.argv[1]

rlist = []
length = 0
with open(fn) as fp:
    for line in fp:
        line = line.strip()
        if line[0] == '>':
            if length != 0:
                rlist.append(length)
                length = 0
            continue
        length += len(line)    

#get the last one
rlist.append(length)

if len(set(rlist)) != 1:
    print('Error! the fasta file show inconstent length! Here is the lengths:')
    print('\n'.join(map(str, rlist)))

else:
    print('fasta length is consistent. length is %s' % rlist[0])
