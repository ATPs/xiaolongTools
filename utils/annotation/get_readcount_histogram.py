#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1],"r")
info = fp.readlines()
fp.close()

counts = {}
ids = []
for line in info:
	words = line.split()
	cov = int(words[2])
	try:
		counts[cov] += 1
	except KeyError:
		ids.append(cov)
		counts[cov] = 1
ids.sort()
for id in ids:
	print str(id) + "\t" + str(counts[id])
