#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1],"r")
info = fp.readlines()
fp.close()

repeats = {}
scafs = []
for line in info:
	words = line.split()
	scafinfo = words[0].split(":")
	scaf = scafinfo[0]
	start = int(scafinfo[1].split("-")[0]) - 1
	end = int(scafinfo[1].split("-")[1])
	if not scaf in scafs:
		scafs.append(scaf)
		repeats[scaf] = [[start,end]]
	else:
		repeats[scaf].append([start,end])

fp = open(sys.argv[2],"r")
info = fp.read().split(">")
fp.close()

scafseqs = {}
for item in info:
	if item:
		lines = item.split("\n")
		name = lines[0]
		seq = ""
		for line in lines[1:]:
			seq += line
		scafseqs[name] = seq

os.system("mkdir potential_repeats")
for scaf in scafs:
    try:
        scafseqs[scaf]
        for repeat in repeats[scaf]:
            if repeat[1] > len(scafseqs[scaf]):
                repeat[1] = len(scafseqs[scaf])
            rp1 = open("potential_repeats/" + scaf + "_" + str(repeat[0]) + "-" + str(repeat[1]), "w")
            rp1.write(">" + scaf + ":" + str(repeat[0]) + "-" + str(repeat[1]) + "\n")
            rp1.write(scafseqs[scaf][repeat[0]:repeat[1]] + "\n")
            rp1.close()
    except:
        pass
