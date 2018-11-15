#!/usr/bin/python
import os, sys, string

# usage: demux.py [index_list] [read_fastq] [index_fastq] [read_num]
fp = open(sys.argv[1],"r")
info = fp.readlines()
fp.close()
libs = {}
seqs = {}
for line in info:
	words = line.split()
	libs[words[1]] = words[0]
	seqs[words[1]] = []
seqs["other"] = []

fp1 = open(sys.argv[2],"r")
fp2 = open(sys.argv[3],"r")

item = ""
index = ""
for count, line1 in enumerate(fp1):
	line2 = fp2.readline()
	if count % 4 == 0:
		if index:
			for barcode in list(libs.keys()):
				match = 0
				for countc in range(len(index)):
					if index[countc] == barcode[countc]:
						match += 1
				if match >= 6:
					new_index = barcode
					break
			else:
				new_index = "other"
			seqs[new_index].append(item)
		item = line1				
	else:
		if count%4 == 1:
			index = line2[:-1]
		item += line1
for id in list(libs.keys()):
	rp = open(libs[id] + "_R" + sys.argv[4] + ".fq", "w")
	for item in seqs[id]:
		rp.write(item)
	rp.close()
rp = open("others_R" + sys.argv[4] + ".fq", "w")
for item in seqs["other"]:
	rp.write(item)
rp.close()
