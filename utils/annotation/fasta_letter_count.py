#!/usr/bin/python
import os, sys, string

fp = open(sys.argv[1],"r")
info = fp.read().split(">")
fp.close()

for item in info:
    if item:
        lines = item.split("\n")
        name = lines[0]
        seq = ""
        for line in lines[1:]:
            seq += line
        count = 0
        for char in seq:
            if char != "N":
                count += 1
        print( name + "\t" + str(len(seq)) + "\t" + str(count))
