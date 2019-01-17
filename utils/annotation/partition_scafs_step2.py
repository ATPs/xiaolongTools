#!/usr/bin/python
import os, sys, string

segments = []
fp = open(sys.argv[1],"r").readlines()
for line in fp:
	words = line.split()
	if int(words[1]) >= 100:
		segments.append(words)	
	else:
		seginfo = words[0].split(":")
		segstart = seginfo[1].split("-")[0]
		segend = seginfo[1].split("-")[1]
		if segments[-1][0].split(":")[0] == seginfo[0]:
			segments[-1][0] = segments[-1][0].split("-")[0] + "-" + str(segend)
			segments[-1][1] = str(int(segments[-1][1]) + int(words[1]))
		else:
			segments.append(words)
for segment in segments:
	print segment[0] + "\t" + segment[1]
