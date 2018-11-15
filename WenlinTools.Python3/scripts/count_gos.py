#!/usr/bin/env python
import os, sys, string
import scipy
from scipy import stats

anns = {}
types = {}
fp = open("/project/biophysics/Nick_lab/wli/sequencing/scripts/data/go.annotation","r")
info = fp.readlines()
fp.close()
for line in info:
    words = line.split("\t")
    anns[words[0]] = words[2]
    types[words[0]] = words[1]

fp = open("/project/biophysics/Nick_lab/wli/sequencing/scripts/data/go.parents","r")
info = fp.readlines()
fp.close()
parents = {}
for line in info:
    words = line.split()
    parents[words[0]] = words[1].split(",")

prots = {}
#usage *.py query.gos background.gos
fp = open(sys.argv[1],"r")
info = fp.readlines()
fp.close()
n = 0
gos = {}
for line in info:
    n += 1
    words = line.split()[1].split(",")
    current = []
    for word in words:
        current.append(word)
        try:
            prots[word].add(line.split()[0])
        except KeyError:
            prots[word] = set([line.split()[0]])
        try:
            current += parents[word]
            for parent in parents[word]:
                try:
                    prots[parent].add(line.split()[0])
                except KeyError:
                    prots[parent] = set([line.split()[0]])
        except KeyError:
            pass
    for go in set(current):
        try:
            gos[go] += 1
        except KeyError:
            gos[go] = 1

fp = open(sys.argv[2],"r")
info = fp.readlines()
fp.close()
allgos = {}
N = 0
for line in info:
    N += 1
    words = line.split()[1].split(",")
    current = []
    for word in words:
        current.append(word)
        try:
            current += parents[word]
        except KeyError:
            pass
    for go in set(current):
        try:
            allgos[go] += 1
        except KeyError:
            allgos[go] = 1

print(len(list(gos.keys())))
for go in list(gos.keys()):
    m = gos[go]
    p = float(allgos[go])/N
    P = scipy.stats.binom_test(m,n,p)
    if P < 0.01:
        try:
            print(go + "\t" + str(P) + "\t" + str(m) + " " + str(n) + " " + str(allgos[go]) + " " + str(N) + "\t" + types[go] + "\t" + anns[go] + "\t" + string.join(prots[go],","))
        except KeyError:
            print(go + "\t" + str(P) + "\t" + str(m) + " " + str(n) + " " + str(allgos[go]) + " " + str(N))
