#!/usr/bin/python
import os, sys, string

fp1 = open(sys.argv[1],"r")
fp2 = open(sys.argv[2],"r")

for count, line1 in enumerate(fp1):
        line2 = fp2.readline()
        if count % 4 == 1 and len(line1) >= 150 and len(line2) >= 150:
                print(line1[10:60] + line2[10:60])
