#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1],"r")
scaf2size = {}
for line in fp:
	words = line.split()
	try:
		scaf2size[words[0]] += 1
	except KeyError:
		scaf2size[words[0]] = 1
fp.close()

rp = open("scaf_sizes","w")
for scaf in scaf2size.keys():
	rp.write(scaf + "\t" + str(scaf2size[scaf]) + "\n")
rp.close()
