#!/usr/bin/python
import os, sys, string

# sys.argv[1] is the basename of the input data
# sys.argv[2] is the coverage cutoff (maximum gap fraction allowed) to remove gappy position

gapfrac = float(sys.argv[2])

fp = open(sys.argv[1],"r")
rp = open(sys.argv[1] + "_gap" + sys.argv[2],"w")
for countl, line in enumerate(fp):
	words = line.split()
	if not countl:
		rp.write(line)
		total = len(words) - 2
	else:
		gapcount = 0
		for countw, word in enumerate(words[2:]):
			if word == "-" or word == "N":
				gapcount += 1
		if gapcount <= total * gapfrac:
			rp.write(line)
fp.close()
rp.close()

