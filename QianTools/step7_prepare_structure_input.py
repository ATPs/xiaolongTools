#!/usr/bin/python
import os, sys, string

pops = []
fp = open(sys.argv[1] + ".ind","r")
for line in fp:
	words = line.split()
	if not words[2] in pops:
		pops.append(words[2])
fp.close()

pop2num = {}
rp = open("population_name_and_population_code","w")
for i, pop in enumerate(pops):
	pop2num[pop] = str(i + 1)
	rp.write(pop + "\t" + str(i+1) + "\n")
rp.close()

fp = open(sys.argv[1] + ".ind","r")
id2info = {}
for line in fp:
	words = line.split()
	id2info[words[0]] = pop2num[words[2]]
	id2info[words[0]] = pop2num[words[2]]
	id2info[words[0]] = pop2num[words[2]]
	id2info[words[0]] = pop2num[words[2]]
	id2info[words[0]] = pop2num[words[2]]
fp.close()

fp = open(sys.argv[1],"r")
rp = open(sys.argv[1] + ".struct","w")
for count, line in enumerate(fp):
	words = line.split()
	if not count:
		sps = words[2:]
		sp2snp1 = {}
		sp2snp2 = {}
		for sp in sps:
			sp2snp1[sp] = []
			sp2snp2[sp] = []
	else:
		nts = []
		for word in words[2:]:
			for subword in word.split(","):
				if subword != "-":
					if not subword in nts:
						nts.append(subword)
		
		ntA = nts[0]
		scaf = words[0]
		for countw, word in enumerate(words[2:]):
			sp = sps[countw]
			nt1 = word.split(",")[0]
			nt2 = word.split(",")[1]
			if nt1 == "-":
				sp2snp1[sp].append("-9")
			elif nt1 == ntA:
				sp2snp1[sp].append(" 0")
			else:
				sp2snp1[sp].append(" 1")
			if nt2 == "-":
				sp2snp2[sp].append("-9")
			elif nt2 == ntA:
				sp2snp2[sp].append(" 0")
			else:
				sp2snp2[sp].append(" 1")
fp.close()

countind = 0
countpos = 0
for sp in sps:
	if not countpos:
		countpos = len(sp2snp1[sp])
	countind += 1
	rp.write(sp.ljust(20) + id2info[sp] + "  0  " + string.join(sp2snp1[sp]," ") + "\n")
	rp.write(sp.ljust(20) + id2info[sp] + "  0  " + string.join(sp2snp2[sp]," ") + "\n")
rp.close()

print "number_of_individuals\t" + str(countind)
print "number_of_positions\t" + str(countpos)
