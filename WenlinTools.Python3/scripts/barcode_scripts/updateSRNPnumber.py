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
import requests as rq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        ID=sys.argv[1]
    except:
        print("Usage: *.py SRNPnumber", file=sys.stderr)
        sys.exit()


    url = 'http://janzen.sas.upenn.edu/Wadults/resultsexpressVOUCHB.lasso'

    json = {'submitButtonName': 'SUBMIT',
            'voucher': ID}

    r = rq.post(url, data=json)

    sp = cmn.find_between(r.content, 'species:<b>', '</b>').strip()

    if '<title>' not in sp:
        print(sp)



