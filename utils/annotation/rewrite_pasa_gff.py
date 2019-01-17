#!/usr/bin/python
import os,sys,string

fp = open(sys.argv[1],"r")
info = fp.readlines()
fp.close()

gff = {}
starts = {}
ends = {}
scafs = {}
dirs = {}
ids = []
for line in info:
    words = line.split('\t')
    if len(words) > 8:
        try:
            transcript_id = words[8].split(';')[1].split()[1].replace('"','')
        except:
            print('error for line', line)
            exit()
        try:
            gff[transcript_id].append(line[:-1])
            starts[transcript_id].append(int(words[3]))
            ends[transcript_id].append(int(words[4]))
        except KeyError:
            gff[transcript_id] = [line[:-1]]
            starts[transcript_id] = [int(words[3])]
            ends[transcript_id] = [int(words[4])]
            scafs[transcript_id] = words[0]
            dirs[transcript_id] = words[6]
            ids.append(transcript_id)
            
for transcript_id in ids:
    print (scafs[transcript_id] + "\tPASA\ttranscript\t" + str(min(starts[transcript_id])) + "\t" + str(max(ends[transcript_id])) + "\t.\t" + dirs[transcript_id] + "\t.\tID=" + transcript_id)
    for line in gff[transcript_id]:
        print (scafs[transcript_id] + "\tPASA\texon\t" + string.join(line.split("\t")[3:],"\t"))
