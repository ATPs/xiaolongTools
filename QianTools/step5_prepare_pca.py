#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1],"r") # the file contains the sample information, see example.list 
id2info = {}
for line in fp:
	words = line.split()
	id2info[words[0]] = [words[1],words[2]]
fp.close()

fp = open(sys.argv[2],"r") # the file that contains the sizes of each scaffold
scafsize = {}
for line in fp:
	words = line.split()
	scafsize[words[0]] = int(words[1])
fp.close()

fp = open(sys.argv[3],"r")
scafs = []
scafset = set([])
for countl, line in enumerate(fp):
	if countl:
		words = line.split()
		if not words[0] in scafset:
			scafs.append(words[0])
			scafset.add(words[0])
fp.close()

before = 0
beforescaf = {}
for scaf in scafs:
	beforescaf[scaf] = before
	before += scafsize[scaf] + 100000

fp = open(sys.argv[3],"r")
rp1 = open(sys.argv[3] + ".snp","w")
rp2 = open(sys.argv[3] + ".geno","w")
rp3 = open(sys.argv[3] + ".ind","w")

for count, line in enumerate(fp):
	words = line.split()
	if not count:
		sps = words[2:]
		for sp in sps:
			rp3.write(sp + "\t" + id2info[sp][0] + "\t" + id2info[sp][1] + "\n")
	else:
		nts = []
		for word in words[2:]:
			for subword in word.split(","):
				if subword != "-":
					if not subword in nts:
						nts.append(subword)
		if len(nts) == 2:
			newline = []
			nt1 = nts[0]
			nt2 = nts[1]
			scaf = words[0]
			num = beforescaf[scaf] + int(words[1])
			rp1.write("snp" + str(count) + "\t1\t0.0\t" + str(num) + "\t" + nt1 + "\t" + nt2 + "\n")
			for word in words[2:]:
				count1 = 0
				count2 = 0
				for subword in word.split(","):
					if subword == nt1:
						count1 += 1
					elif subword == nt2:
						count2 += 1
				if count1 + count2 < 2:
					newline.append("9")
				else:
					newline.append(str(count1))
			rp2.write(string.join(newline,"") + "\n")
rp1.close()
rp2.close()
rp3.close()
