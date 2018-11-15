import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


result_file=sys.argv[1]


dic={}
dic['A']='T'
dic['T']='A'
dic['C']='G'
dic['G']='C'
dic['N']='N'
dic['-']='N'

lines = cmn.file2lines(result_file)

if len(lines) == 0:
    print('no blast result!')
    sys.exit()


sorted_result = sorted(lines, key = lambda x: int(x.split()[4]), reverse=True)
bestline = sorted_result[0]

#the result included 20 primers, so need to check this
qseqid, sseqid, evalue, pident, length, qlen, qstart, qend, slen, sstart, send, qseq, sseq = bestline.strip().split()

sstart = int(sstart)
send = int(send)

isReverse = False
if sstart > send:
    sstart, send = send, sstart
    isReverse = True

if sstart <= 20 and send >= 678:
    isComplete = True
else:
    isComplete = False


barcode = qseq
if isReverse:
    barcode = ''.join([dic[char] for char in barcode[::-1]])


if sstart <= 20:
    leftShift = 20 - sstart + 1
    barcode = barcode[leftShift:]
else:
    add_Ns = sstart - 20
    barcode = '%s%s' % ('N'*add_Ns, barcode)

if send >= 678:
    rightShift = send - 678
    barcode = barcode[:-rightShift]
else:
    add_Ns = 678 - send
    barcode = '%s%s' % (barcode, 'N'*add_Ns)


dn = 'denovo_barcode.fa'
fasta = '>denovo_barcode\n%s\n' % barcode
cmn.write_file(fasta, dn)
