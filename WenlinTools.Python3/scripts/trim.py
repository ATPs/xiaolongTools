#!/usr/bin/python
import os, sys, string

seq = ""
name = ""
for count, line in enumerate(sys.stdin):
	if count % 4 == 0:
		name = line[:-1]
	elif count % 4 == 1:
		seq = line
		total = len(seq)
	elif count % 4 == 3:
		start = 0
		end = 0
		for char in line:
			if char > '5':
				break
			start += 1
		for char in line[::-1]:
			if char > '5':
				break
			end -= 1
		if total + end - start >= 30:
			print(name)
			print(seq[start:end]) 
			print("+")
			print(line[start:end])	
