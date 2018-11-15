#functions that access uniprot eUtils
import cmn,sys

biopath='/home/wenlin/local/biopython-1.58'

if biopath not in sys.path:
    sys.path.insert(0,biopath)
from Bio import ExPASy
from Bio import SwissProt

import urllib.request, urllib.parse, urllib.error,urllib.request,urllib.error,urllib.parse
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#used by IDmapping
#go to server to get the uniprot ID, return a dict like dict[fromgi]=[togi1,togi2]
def gi2server(fromdb,todb,id_list):
    query='  '.join(id_list)
    content=100*len(query)
    #ID=[]
    url = 'http://www.uniprot.org/mapping/'

    params = {
        'from':fromdb,
        'to':todb,
        'format':'tab',
        'query':query
        }


    data = urllib.parse.urlencode(params)
    request = urllib.request.Request(url, data)
    response = urllib.request.urlopen(request)
    page = response.read(content)

    page=page.split('\n')
    page=page[1:]
    #page like [fromgi\ttogi]
    #adict={}
    #for line in page:
     #   if line!='':
      #      item=line.split()
       #     ID=[item[-1]]
    return page
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the result_list like ['gi\tuniID']
#the mapping tag can be found in http://www.uniprot.org/faq/28
#mapping id between databases
class IDmapping(object):
    def __init__(self,idlist,fromdb='P_GI',todb='ACC',bundle=200):
        if type(idlist) is not list:#single gi
            idlist=[idlist]
        self.fromid=idlist#the same order as input
        self.mapped=self.mapping(fromdb,todb,bundle)
        self.toid=self.gettoid()


    def mapping(self,fromdb,todb,bundle):
        id_list=self.fromid
        if len(id_list)<bundle:
            result_list=gi2server(fromdb,todb,id_list)
        else:
            result_list=[]
            increase=bundle
            total_num=len(id_list)
            times=total_num/increase
            #go to the server
            for i in range(times+1):
                temp_list=id_list[i*increase:(i+1)*increase]
                result_list+=gi2server(fromdb,todb,temp_list)
        result_list=cmn.del_empty(result_list)
        #result list like ['gi\tuniID']
        adict={}
        for line in result_list:
            item=line.split()
            key=item[0]
            value=item[1]
            if key not in adict:
                adict[key]=[value]
            else:#key exist
                adict[key]+=[value]

        return adict

    def gettoid(self):
        adict=self.mapped
        idlist=[]
        for i in adict:
            temp=adict[i]#each mapped record is a list
            for v in temp:
                idlist+=[v]
        return idlist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class IDsummary(object):

    def __init__(self,idlist):
        self.fromid=idlist
        self.esum=self.getesum()

    def getesum(self):
        id_list=self.fromid
        record=[]
        record+=self.eachget(id_list)
        return record

    def eachget(self,id_list):
        a=[]
        for eachid in id_list:
            try:
                record=SwissProt.read(ExPASy.get_sprot_raw(eachid))
                #print 'testing\n'
                a+=[record]
            except:#something to do if no summary found
                print('something wrong with this id:%s\n' % eachid)#here do nothing
        return a








