#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1],"r")
info = fp.read().split(">")
fp.close()
lens = []
totallen = 0
rp = open(sys.argv[2],"w")
for item in info:
	lines = item.split("\n")
	name = lines[0]
	seq = ""
	for line in lines[1:]:
		seq += line
	lens.append(len(seq))
	totallen += len(seq)
	rp.write(name + "\t" + str(len(seq)) + "\n")
rp.close()

lens.sort()
current = 0
for length in lens:
	current += length
	if current >= totallen * 0.5:
		print(length)
		break
