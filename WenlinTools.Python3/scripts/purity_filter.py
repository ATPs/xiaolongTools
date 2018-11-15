#!/usr/bin/python
import os, sys, string

fn = sys.argv[1]
fp = open(sys.argv[1],"r")
label = fn.split('/')[-1].split('.')[0].split('_')[-1][-1]

for count, line in enumerate(fp):
	if count % 4 == 0:
		words = line.split()
		subwords1 = words[0].split(":")
		subwords2 = words[1].split(":")
		if subwords2[1] == "N":
			get = 1
			#print "@" + string.join(subwords1[4:7],":") + "/" + label
			print('%s/%s' % (words[0], label))
		else:
			get = 0
	elif get:
		print(line[:-1])
