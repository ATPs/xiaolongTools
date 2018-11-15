import sys
fn=sys.argv[1]
cutoff_ratio=float(sys.argv[2])
dic={}
f=open(fn)
count=0
for i in f:
    if '>' not in i:
        count=count+1
        for k in range(len(i)-1):
            if i[k]!='-' and i[k]!='N':
                try:
                    dic[k]=dic[k]+1
                except:
                    dic[k]=1
f.close()
pos=[key for key in dic.keys() if dic[key]>=cutoff_ratio*count]
fw=open(fn.split('.')[0]+'_rmGap_'+str(cutoff_ratio)+'.fasta','w')
f=open(fn)
for i in f:
    if '>' in i:
        fw.write(i)
    else:
        fw.write(''.join([i[k] for k in pos]))
        fw.write('\n')
f.close()
fw.close()

