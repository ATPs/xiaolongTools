#!/usr/bin/python
import sys

fp = open(sys.argv[1],"r")
info = fp.readlines()
fp.close()

counts = {}
ids = []
for line in info:
    words = line.split()
    cov = int(round(float(words[4]),0))
    try:
        counts[cov] += int(words[2])
    except KeyError:
        ids.append(cov)
        counts[cov] = int(words[2])
ids.sort()
for i in ids:
    print (str(i) + "\t" + str(counts[i]))