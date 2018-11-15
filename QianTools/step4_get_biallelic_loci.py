#!/usr/bin/python
import os, sys, string

# CUTOFF1 is how many times each allele is required to show up
# CUTOFF2 is how many specimens each allele is required to show up
# these cutoffs need to change depending on the expected size of each major populations. If there are about 10 samples per population, the below cutoff may be good. If there are about 5 samples per population, you may use 4 and 3. But CUTOFF1 should always be larger than 3, and CUTOFF2 should be always larger than 2. 
CUTOFF1 = 8
CUTOFF2 = 6

fp = open(sys.argv[1], "r")
rp = open(sys.argv[1] + "_biallelic", "w")
for count, line in enumerate(fp):
	words = line.split()
	if not count:
		sps = []
		for word in words[2:]:
			sp = word[:-2]
			sps.append(sp)
		rp.write(words[0] + "\t" + words[1] + "\t" + string.join(sps,"\t") + "\n")
	else:
		chars = set([])
		char2count = {}
		char2sp = {}
		sp2chars = {}
		for countw, word in enumerate(words[2:]):
			sp = sps[countw]
			if word == "N":
				nt = "-"
			else:
				nt = word
			try:
				sp2chars[sp].append(nt)
			except KeyError:
				sp2chars[sp] = [nt]
			if nt != "-":
				chars.add(nt)
				try:
					char2count[nt] += 1
				except KeyError:
					char2count[nt] = 1
				try:
					char2sp[nt].add(sp)
				except KeyError:
					char2sp[nt] = set([sp])
		if len(chars) == 2:
			min1 = 1000
			min2 = 1000
			for char in char2count.keys():
				if char2count[char] < min1:
					min1 = char2count[char]
			for char in char2sp.keys():
				if len(char2sp[char]) < min2:
					min2 = len(char2sp[char])
			if min1 >= CUTOFF1 or min2 >= CUTOFF2:
				newline = [words[0],words[1]]
				for sp in sps:
					newline.append(string.join(sp2chars[sp],","))
				rp.write(string.join(newline,"\t") + "\n")
fp.close()
rp.close()
