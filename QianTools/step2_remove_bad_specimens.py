#!/usr/bin/python
import os, sys, string

# sys.argv[1] is the basename of the input data
# sys.argv[2] is the coverage cutoff to remove bad (incomplete) samples

fp = open(sys.argv[1] + ".sample_coverage","r")
remove = set([])
for line in fp:
	words = line.split()
	if float(words[3]) < float(sys.argv[2]):
		remove.add(words[0])
fp.close()

fp = open(sys.argv[1] + ".list","r")
sps = []
newsps = []
for line in fp:
	words = line.split()
	sps.append(words[0])
	if not words[0] in remove:
		newsps.append(words[0] + "_1")
		newsps.append(words[0] + "_2")
fp.close()

fp = open(sys.argv[1] + ".data","r")
rp = open(sys.argv[1] + "_nobad" + sys.argv[2],"w")
rp.write("scaffold\tposition\t" + string.join(newsps,"\t") + "\n")
for line in fp:
	words = line.split()
	newwords = words[:2]
	for countw, word in enumerate(words[2:]):
		sp = sps[countw/2]
		if not sp in remove:
			newwords.append(word)
	rp.write(string.join(newwords,"\t") + "\n")
fp.close()
rp.close()

