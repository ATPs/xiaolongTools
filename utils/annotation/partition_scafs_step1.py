#!/usr/bin/python
import os, sys, string

scfs = []
fp = open(sys.argv[1],"r").readlines()
for line in fp:
	words = line.split()
	scaf = words[0]
	scfs.append(scaf)

scafseqs = {}
fp = open(sys.argv[2],"r").read().split(">")
for item in fp:
	if item:
		lines = item.split("\n")
		scaf = lines[0]
		seq = ""
		for line in lines[1:]:
			seq += line
		scafseqs[scaf] = seq

for scf in scfs:
	subseqs = scafseqs[scf].split("N")
	start = 0
	for subseq in subseqs:
		start += 1
		end = start + len(subseq)
		if end == start:
			continue
		elif end - start <= 100:
			print scf + ":" + str(start) + "-" + str(end) + "\t" + str(end - start)
		else:
			substart = 0
			for i in range(len(subseq)/100):
		                substart = i*100
				print scf + ":" + str(start + substart) + "-" + str(start + substart + 100) + "\t100"
			print scf + ":" + str(start + substart + 100) + "-" + str(end) + "\t" + str(end - start - substart - 100)
		start += len(subseq)
