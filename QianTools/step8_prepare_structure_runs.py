#!/usr/bin/python
import os, sys, string

name = sys.argv[1]
jobhead = "#!/bin/bash\n#SBATCH -p long\n#SBATCH -n 1\n#SBATCH -N 1\n#SBATCH --mem=8g\n"

cwd = os.getcwd()
os.system("mkdir run")
for size in ["2","3","4","5","6"]:
	os.system("mkdir run/" + size)
	os.system("cp extraparams run/" + size)
	countp = os.popen("head -1 " + name + ".struct | wc")
	countinfo = countp.readline()
	countwords = countinfo.split()
	count = int(countwords[1]) - 3
	fp = open("mainparams","r")
	rp = open("run/" + size + "/mainparams","w")
	for line in fp:
		words = line.split()
		if len(words) > 2:
			if words[0] == "#define" and words[1] == "NUMLOCI":
				rp.write(words[0] + " " + words[1] + "    " + str(count) + " " + string.join(words[3:]," ") + "\n")
			else:
				rp.write(line)
		else:
			rp.write(line)
	fp.close()
	rp.close()
	
	os.system("cp " + name + ".struct run/" + size)
	rp = open(name + "_" + size + ".job","w")
	rp.write(jobhead)
	rp.write("#SBATCH -o " + name + "_" + size + ".log\n\n")
	rp.write("cd " + cwd + "run/" + size + "\n")
	rp.write("/projects/omics/hesperia/STRUCTURE/console/structure -K " + size + " -i " + name + ".struct -o output_" + size + "\n")
	rp.close()
