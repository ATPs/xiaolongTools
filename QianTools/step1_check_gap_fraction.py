#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1] + ".list","r")
sps = []
sp2count = {}
for line in fp:
	words = line.split()
	sps.append(words[0])
	sp2count[words[0]] = 0
fp.close()

total = 0
fp = open(sys.argv[1] + ".data","r")
for line in fp:
	words = line.split()
	total += 1
	for countw, word in enumerate(words[2:]):
		if countw % 2 == 0:
			sp = sps[countw/2]
			if word != "-" and word != "N":
				sp2count[sp] += 1
fp.close()

rp = open(sys.argv[1] + ".sample_coverage","w")
for sp in sps:
	rp.write(sp + "\t" + str(total) + "\t" + str(sp2count[sp]) + "\t" + str(float(sp2count[sp])/total) + "\n")
rp.close()
