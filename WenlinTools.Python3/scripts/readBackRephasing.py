#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_options():
    import optparse, os
    parser = optparse.OptionParser()

    #help info
    usage= '\nrunhhblits.py -i query.fasta [-o result.hhr] [-d database]\n\n'
    parser = optparse.OptionParser(usage=usage)


    #define your option here
    parser.add_option('-i', dest='infile', help='input sequence file (required)')

    options,args = parser.parse_args()

    #checking possible typo
    #no args
    if len(sys.argv) == 1 and len(args) == 0:
        os.system(sys.argv[0] + ' -h')
        sys.exit()
    #args format problem
    if len(args) != 0:
        print("error: incorrect number of arguments.\n")
        print("please use \'seq2ref.py -h\' to print the help message.\n")
        sys.exit()
    if not options.infile:
        print("error: No sequence file input!\n")
        sys.exit()

    return options



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
def main():
    True

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()






