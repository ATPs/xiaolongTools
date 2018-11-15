#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#main
def main():
    True

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=='__main__':
    import sys
    alist=sys.stdin.readlines()
    #b should contain %s
    b=sys.argv[1]

    for a in alist:
        print(b.replace('%s', a.strip()))

