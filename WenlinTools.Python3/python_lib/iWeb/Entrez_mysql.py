import sys

biopath='/home/wenlin/local/biopython-1.58'

if biopath not in sys.path:
    sys.path.insert(0,biopath)
from database.mysql2 import Entrez as myEn
import _mysql_exceptions

from Bio import Entrez
Entrez.email='wenlin.work@hotmail.com'

import time
from datetime import datetime

import cmn

from Bio.Entrez.Parser import NotXMLError

from iObject.shelves import shelveIO as sh


def ifconnect(testgi='254780193'):
    """
    Used to test whether NCBI block the ip
    return True or False. if good, return True
    if it is good, pause 5 secs.
    """
    h=Entrez.elink(dbfrom='protein',db='pubmed',id=testgi)
    try:
        Entrez.read(h)
        time.sleep(5)
        return True
    except NotXMLError:
        return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def map_dict(fromdb,todb,adict,record,linkname=None, recording=None):
    if linkname:
        dbname = "Entrez_%s" % linkname
    else:
        dbname = "Entrez_%s2%s" % (fromdb, todb)
    refdb = myEn(dbname)
    ct = cmn.today()
    #deal with record
    tmpdict={}
    for entry in record:
        key=entry["IdList"][0]
        #print 'dealing with id: '+key
        value=[]
        for linked_db in entry["LinkSetDb"]:#linked_db like [DbTo, Link, LinkName]
            for linked_id in linked_db["Link"]:
                value+=[linked_id["Id"]]
        tmpdict[key]=value
        refdb.insert([key, ','.join(value), ct])
    #if recording:
    #    save_local(tmpdict,"%s_%s" % (fromdb,todb))
    adict.update(tmpdict)
    refdb.close()
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def detail_map_dict(adict,ddict,record):
    #deal with record
    for entry in record:
        key1=entry["IdList"][0]
        #print 'dealing with id: '+key
        value={}#used in ddict
        a_value=[]#used in adict
        for linked_db in entry["LinkSetDb"]:#linked_db like [DbTo, Link, LinkName]
            key2=linked_db["LinkName"]
            dlist=[]
            if True:
                for linked_id in linked_db["Link"]:
                    dlist+=[linked_id["Id"]]
                    a_value+=[linked_id["Id"]]
                value[key2]=dlist
        ddict[key1]=value
        adict[key1]=a_value
    return adict,ddict
        # ddcit looks like adict[gi]={protein_pubmed:***;protein_pubmed_refseq}
        # adict just adict[gi]=[pubid1,pubid2]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#exclusion used in the 'LinkName'
def detail_goserver(r_dict,d_dict,fromdb,todb,alist):
    try:#the normal way
        record=Entrez.read(Entrez.elink(dbfrom=fromdb,db=todb,id=alist))
        r_dict,d_dict=detail_map_dict(r_dict,d_dict,record)
        del record
        #print 'I am normal!\n'
    except:
        #first try do it one by one
        for aa in alist:
            try:
                print('try it one by one, here I am in %s .\n' % aa)
                record=Entrez.read(Entrez.elink(dbfrom=fromdb,db=todb,id=aa))
                r_dict,d_dict=detail_map_dict(r_dict,d_dict,record)
                del record
                #print 'I succeed at my first trial!\n'
            except:
                #if one still fail, try sometimes
                print('I need try more times on %s .\n' % aa)
                flag=False
                count=0#times trying
                while(flag==False and count<3):
                    count+=1
                    print('trying my %s times...\n' % count)
                    try:
                        time.sleep(2)#interval between times
                        record=Entrez.read(Entrez.elink(dbfrom=fromdb,db=todb,id=aa))
                        r_dict,d_dict=detail_map_dict(r_dict,d_dict,record)
                        del record
                        flag=True
                        #print 'yes, I have done it!\n'
                    except:
                        flag=False
                #well, the server doesn't like me...
                if flag==False:
                    print('something wrong with this id: %s\n' % aa)

    for key in r_dict:
        r_dict[key]=cmn.unify(r_dict[key])
    return r_dict,d_dict


#use the shelve database to do the fetch
#if record older than duration days, do not retrieve
#the dict is designed to be d[fromid]={"date":2013-03-14,"ids":["23232"]}
#TODO: need to debug this
def local_fetch(r_dict,idlist,tableName,duration=10):
    adict = {}
    my = Entrez(tableName)
    ctime=datetime.now()
    bundle = 100
    for i in range(len(idlist)/bundle + 1):
        tmp = idlist[i*bundle: (i+1)*bundle]
        records = my.select("where fromid in (%s)", ','.join(tmp))
        for record in records:
            fromid, toid, date = record
            #TODO: debug the date
            if (ctime - date).days < duration:
                adict[fromid] = toid.split(",")

    remainingIDlist=[]
    for ID in idlist:
        if ID not in adict:
            remainingIDlist+=[ID]
    return adict,remainingIDlist
#thedict is like dict[fromid]=[toidlist]
#final dict is designed to be d[fromid]={"date":2013-03-14,"ids":["23232"]}
def save_local(thedict,filename):
    base='/data/wenlin_database/shelves/Entrez/'
    ctime=datetime.now()
    newdict={}
    for ID in thedict:
        newdict[ID]={
                "date": "%s-%s-%s" % (ctime.year,ctime.month,ctime.day),
                "ids": thedict[ID]
                }
    my=sh(base+filename)
    my.update(newdict)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TODO: add something related to local fetch here
def visit_Entrez(fromdb="protein", todb="pubmed", alist=[], linkname=None):
    record = []
    count = 0
    while (record == [] and count < 3):
        if linkname:
            record=Entrez.read(Entrez.elink(dbfrom=fromdb,db=todb,id=alist,linkname=linkname))
        else:
            record=Entrez.read(Entrez.elink(dbfrom=fromdb,db=todb,id=alist))
        if record == []:
            count += 1
            print("Weird info from Entrez server, try again (%s trial)." % count)
            time.sleep(2)
    if record == []:
        print("Error in Entrez server communication for:\n%s" % ",".join(alist), file=sys.stderr)

    return record

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def goserver(r_dict,fromdb,todb,alist,linkname=None, recr=False, duration=15, recording=False,times=0):
    if len(alist) == 0:
        return {}
    ct = datetime.now()
    #first check local mysql database
    if linkname:
        dbname = "Entrez_%s" % linkname
    else:
        dbname = "Entrez_%s2%s" % (fromdb, todb)
    refdb = myEn(dbname)
    try:
        refdb.createTable(['%s INT UNSIGNED, %s TEXT, update_date DATE, PRIMARY KEY (%s)' % (fromdb, todb, fromdb)])
    except _mysql_exceptions.OperationalError:
        pass

    for a in alist:
        mapping = refdb.select("WHERE %s=%s" % (fromdb, a))
        if len(mapping) != 0:
            date = mapping[0][-1]
            if (ct-date).days < duration:
                r_dict[a] = mapping[0][1].split(",")

    remaining = set(alist) - set(r_dict.keys())
    try:#if everything is good
        record = visit_Entrez(fromdb=fromdb, todb=todb, alist=remaining, linkname=linkname)
        r_dict = map_dict(fromdb,todb,r_dict,record,linkname=linkname,recording=recording)
    except:
        if time == 1:#take everyone and fetch
            for i in remaining:
                record = visit_Entrez(fromdb=fromdb, todb=todb, alist=[i], linkname=linkname)
                r_dict = map_dict(fromdb,todb,r_dict,record,linkname=linkname,recording=recording)
        else:
            times += 1
            midpoint=len(alist)/2
            aalist=alist[:midpoint]
            ablist=alist[midpoint:]
            r_dict=goserver(r_dict,fromdb,todb,aalist,linkname=linkname,recr=True,recording=recording, times=times)
            r_dict=goserver(r_dict,fromdb,todb,ablist,linkname=linkname,recr=True,recording=recording, times=times)

    return r_dict




    refdb.close()
#exclusion used in the 'LinkName'
def old_goserver(r_dict,fromdb,todb,alist,linkname=None, recr=False,recording=False,times=0):
    try:#the normal way
    #if True:
        if recording:
            if recr==False:#not in recursive loop
                r_dict,alist=local_fetch(r_dict,alist,"%s_%s" % (fromdb,todb))
        if alist != []:
            record = visit_Entrez(fromdb=fromdb, todb=todb, alist=alist, linkname=linkname)
        #print record
        r_dict=map_dict(fromdb,todb,r_dict,record,recording=recording)
        del record
        #print 'I am normal!\n'
    #else:
    except:
        #first try do it one by one
        #for aa in alist:
        if len(alist)==1:
            aa=alist[0]
            #try:
            #    print 'try it one by one, here I am in %s .\n' % aa
            #    record=Entrez.read(Entrez.elink(dbfrom=fromdb,db=todb,id=aa))
            #    r_dict=map_dict(fromdb,todb,r_dict,record,exclusion)
            #    del record
                #print 'I succeed at my first trial!\n'
            #except:
                #if one still fail, try sometimes
            if True:
                print('I need try more times on %s .\n' % aa)
                flag=False
                count=0#times trying
                while(flag==False and count<3):
                    count+=1
                    print('trying my %s times...\n' % count)
                    try:
                        time.sleep(3)#interval between times
                        record = visit_Entrez(fromdb=fromdb, todb=todb, alist=alist, linkname=linkname)

                        r_dict=map_dict(fromdb,todb,r_dict,record,recording=recording)
                        del record
                        flag=True
                        #print 'yes, I have done it!\n'
                    except:
                        flag=False
                #well, the server doesn't like me...
                if flag==False:
                    print('something wrong with this id: %s\n' % aa)
        elif times < 3:
            #if the length is not one yet, divide and do again
            #only go into recursion for 3 loop
            times += 1
            try:
                midpoint=len(alist)/2
                aalist=alist[:midpoint]
                ablist=alist[midpoint:]
                r_dict=goserver(r_dict,fromdb,todb,aalist,linkname=linkname,recr=True,recording=recording, times=times)
                r_dict=goserver(r_dict,fromdb,todb,ablist,linkname=linkname,recr=True,recording=recording, times=times)
            except RuntimeError as theError:
                if theError == "maximum recursion depth exceeded":
                    for item in alist:
                        r_dict=goserver(r_dict,fromdb,todb,[item],linkname=linkname,recr=True,recording=recording)
        else:#more than 3 times
            for item in alist:
                r_dict=goserver(r_dict,fromdb,todb,[item],linkname=linkname,recr=True,recording=recording, times=4)


    for key in r_dict:
        r_dict[key]=cmn.unify(r_dict[key])
    return r_dict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#use exclusion to exclude the database linked
class IDmapping(object):
    def __init__(self,idlist,fromdb="protein",todb="pubmed",bundle=2000, linkname=None, recording=False):
        self.fromid=idlist
        self.mapped=self.mapping(fromdb,todb,bundle,recording,linkname)
        self.toid=self.gettoid()

    def mapping(self,fromdb,todb,bundle,recording,linkname):
        id_list=self.fromid
        result_dict={}
        if len(id_list)<bundle:
            result_dict=goserver(result_dict,fromdb,todb,id_list,linkname=linkname,recording=recording)
        else:
            increase=bundle
            total_num=len(id_list)
            times=total_num/increase
            #go to the server by Biopython
            for i in range(times+1):
                temp_list=id_list[i*increase:(i+1)*increase]
                result_dict=goserver(result_dict,fromdb,todb,temp_list,linkname=linkname, recording=recording)
        return result_dict

    def gettoid(self):
        adict=self.mapped
        idlist=[]
        for i in adict:
            temp=adict[i]#each mapped record is a list
            idlist+=temp
            #for v in temp:
            #    idlist+=[v]
        idlist=cmn.unify(idlist)
        return idlist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class IDsummary(object):
    def __init__(self,idlist,fromdb,bundle=200):
        self.fromid=idlist
        self.esum=self.getesum(fromdb,bundle)

    def getesum(self,checkdb,bundle):
        id_list=self.fromid
        if len(id_list)<bundle:
            if False:
                print(id_list)
                record=Entrez.read(Entrez.esummary(db=checkdb,id=','.join(id_list)))
                print(record)
            else:
                record=[]
                record+=self.eachget(id_list,checkdb)
        else:
            record=[]
            increase=bundle
            total_num=len(id_list)
            times=total_num/increase
            #go to the server by Biopython
            for i in range(times+1):
                temp_list=id_list[i*increase:(i+1)*increase]
                try:
                    record+=Entrez.read(Entrez.esummary(db=checkdb,id=','.join(temp_list)))
                except:
                    record+=self.eachget(temp_list,checkdb)
        return record

    def eachget(self,id_list,checkdb):
        a=[]
        for eachid in id_list:
            count = 0
            record = []
            while(count<3 and record==[]):
                try:
                    record=Entrez.read(Entrez.esummary(db=checkdb,id=eachid))
                    a+=record
                except:#something to do if no summary found
                    True#here do nothing
                count += 1
                if record == []:
                    print("server contact error, wait for 2 second.. (%s trial)" % count)
                    time.sleep(2)
        if a == []:
            print("Weird communication with Entrez server for:\n%s" % ",".join(id_list), file=sys.stderr)
        return a

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#use exclusion to exclude the database linked
#keep the database linkage label
class detail_IDmapping(object):
    def __init__(self,idlist,fromdb,todb,bundle=200):
        self.fromid=idlist
        #self.mapped=self.mapping(fromdb,todb,bundle,exclusion)
        [self.mapped,self.detail_mapped]=self.detail_mapping(fromdb,todb,bundle)
        self.toid=self.gettoid()

    def detail_mapping(self,fromdb,todb,bundle):
        id_list=self.fromid
        result_dict={}
        detail_dict={}
        if len(id_list)<bundle:
            result_dict,detail_dict=detail_goserver(result_dict,detail_dict,fromdb,todb,id_list)
        else:
            increase=bundle
            total_num=len(id_list)
            times=total_num/increase
            #go to the server by Biopython
            for i in range(times+1):
                temp_list=id_list[i*increase:(i+1)*increase]
                result_dict,detail_dict=detail_goserver(result_dict,detail_dict,fromdb,todb,temp_list)
        return [result_dict,detail_dict] #detail_dict: a dict of dict


    def mapping(self,fromdb,todb,bundle):
        id_list=self.fromid
        result_dict={}
        if len(id_list)<bundle:
            result_dict=goserver(result_dict,fromdb,todb,id_list)
        else:
            increase=bundle
            total_num=len(id_list)
            times=total_num/increase
            #go to the server by Biopython
            for i in range(times+1):
                temp_list=id_list[i*increase:(i+1)*increase]
                result_dict=goserver(result_dict,fromdb,todb,temp_list)
        return result_dict

    def gettoid(self):
        adict=self.mapped
        idlist=[]
        for i in adict:
            temp=adict[i]#each mapped record is a list
            idlist+=temp
            #for v in temp:
            #    idlist+=[v]
        idlist=cmn.unify(idlist)
        return idlist


#only accept one once
class Esearch(object):
    def __init__(self,db='',term=''):
        self.record=self.go_search(db,term)
        self.count=self.record['Count']
        self.ids=self.record['IdList']
    def go_search(self,db,term):
        return Entrez.read(Entrez.esearch(db=db,term=term))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#used by mypost
def epost(fromdb="",idlist=[],WebEnv=None):
    if WebEnv==None:#create new
        result=Entrez.read(Entrez.epost(db=fromdb,id=','.join(idlist)))
    else:#append to existing session
        result=Entrez.read(Entrez.epost(db=fromdb,id=','.join(idlist),WebEnv=WebEnv))
    return result["QueryKey"],result["WebEnv"]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class mypost(object):
    """
    use epost to post id list and use efetch to get back the list
    for all the available actions, please refer to:
        http://www.ncbi.nlm.nih.gov/sites/batchentrez
    for accepted input ID, please refer to:
        http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.chapter2_table1/?report=objectonly
    for available output actions, please refer to (most useful one):
        http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly

    3 output parameters, QueryKey, WebEnv and info

    default parameters are assigned according to gi2fasta function
    """
    def __init__(self, fromdb='protein',idlist=[],todb='protein',rettype='fasta',retmode='text',bundle=2222):

        ######post############
        #post the first 2000
        QueryKey,WebEnv=epost(fromdb,idlist[:bundle])
        self.info=Entrez.efetch(db=todb, query_key=QueryKey, WebEnv=WebEnv, rettype=rettype, retmode=retmode).read()
        #if more than 2000, post the rest by bundle
        k=len(idlist)/bundle
        if k>0:#if more than 2000
            for i in range(k):
                ii=i+1
                alist=idlist[ii*bundle:(ii+1)*bundle]
                QueryKey,WebEnv=epost(fromdb,alist,WebEnv=WebEnv)
                self.info+=Entrez.efetch(db=todb, query_key=QueryKey, WebEnv=WebEnv, rettype=rettype, retmode=retmode).read()
        self.QueryKey=QueryKey
        self.WebEnv=WebEnv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


