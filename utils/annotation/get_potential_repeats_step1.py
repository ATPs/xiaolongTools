#!/usr/bin/python
import os, sys, string

fp = open("segment_coverage_100","r")
info = fp.readlines()
fp.close()

repeats = []
for line in info:
	words = line.split()
	if float(words[3]) > 630:
		if repeats:
			if repeats[-1][0].split(":")[0] == words[0].split(":")[0] and repeats[-1][0].split("-")[1] == words[0].split(":")[1].split("-")[0]:
				repeats[-1][0] = repeats[-1][0].split("-")[0] + "-" + words[0].split("-")[1]
				repeats[-1][2] = (repeats[-1][1] * repeats[-1][2] + int(words[1]) * float(words[3]))/(repeats[-1][1] + int(words[1]))
				repeats[-1][1] = repeats[-1][1] + int(words[1])
			else:
				repeats.append([words[0],int(words[1]),float(words[3])])
		else:
			repeats.append([words[0],int(words[1]),float(words[3])])
for repeat in repeats:
	print repeat[0] + "\t" + str(repeat[1]) + "\t" + str(round(repeat[2],0))
