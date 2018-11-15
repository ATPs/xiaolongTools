#function and class commonly used by other programs
import os
import sys
import time
import pickle
import csv


def today():
    from datetime import datetime as dt
    ct = dt.now()
    return "%s-%s-%s" % (ct.year, ct.month, ct.day)

def html2text(info):
    #from bs4 import BeautifulSoup as bs
    if info.strip().startswith('<html'):
        a = info
    else:
        a = "<html>%s</html>" % info
    return a.text

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#delete the number i element in line
def str_del(line,i):
    a=line[:i]+line[i+1:]
    return a

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#get the info between two labels
#input: label, file pointer
#output: the segment, the moved file pointer
def get_seg(label, f):
    line=f.readline()
    if line=='':
        return '',f
    else:
        while(line!='' and (label not in line)):#go to the 1st label
            line=f.readline()
        if line!='':#if the file still has contents
            info=line#first line
            chk=f.readline()
            while(chk!='' and (label not in chk)):#while not reach the next
                info+=chk
                chk=f.readline()
            #reach new one, move back f pointer
            if chk!='':
                back=-(len(chk)+1)
                f.seek=(back,1)
            #if reach the end, do not move back
        else:#no contents found
            return '',f
    return info,f
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#move file pointer back to the beginning of 'line'
#usually, 'line' is just one line
#and 'line' should be the line which f just read
def back_line(f,line):
    backward=-(len(line)+1)
    f.seek(backward,1)
    return f

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#read the file and put the context into a string
#input:file name
#ouput: a string
#string is easier for replace
def txt_read(name):
    f=open(name)
    info=f.read()
    #info=''
    #for line in f:
    #    info+=line
    f.close()
    return info

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#remove duplication in a list
#input and output:list
def unify(alist):
    try:
        sets=set(alist)#sets have no duplicated element
        newlist=list(sets)
    except:
        newlist=[]
        for element in alist:
            if element not in newlist:
                newlist+=[element]
    return newlist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#delete empty item from the list of string
def del_empty(alist):
    has = True
    while (has):
        try:
            alist.remove('')
        except ValueError:
            return alist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#get the info of the link, save to a parameter
#flag: f or info
def link2info(link,wait=None,flag=None,proxy=True):
    if proxy:
        cmd = "export http_proxy='http://proxy.swmed.edu:3128'; "
    else:
        cmd = ""
    cmd += 'wget -qO- \"'+link+'\"'
    if wait!=None:
        time.sleep(wait)
    f=os.popen(cmd)
    if flag!=None:
        if flag=='f':#only need the pointer
            return f
    lines=f.readlines()
    info=''.join(lines)
    f.close()
    return info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#similar function as reading process in gi2fasta.py
#read a gi list, regardless of filename or real gi. return a list
def getid(argv_list):
    if type(argv_list) is str:
        argv_list=[argv_list]

    id_list=[]
    for ids in argv_list:
        if os.path.isfile(ids):#if it is a file
            info=txt_read(ids)
            info=info.split('\n')
            info=del_empty(info)
            id_list+=info
        else:
            id_list+=[ids]
    return id_list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filesize(name):
    statinfo=os.stat(name)
    size=statinfo.st_size
    return size

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if file not exist or file size equal to 0, return false
#wait for 1s and check whether the size is changed
#if it has file, return True
def filexist(name):
    if not os.path.isfile(name):#if not a file
        return False
    else:
        statinfo=os.stat(name)
        size=statinfo.st_size
        if int(size)==0:#no contents
            return False
        else:
            return True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######
def writting(info,d):
    try:
        info=info.encode('UTF-8')
    except:
        pass

    try:
        d.write(info)
    except UnicodeEncodeError:
        try:
            #try utf-8 coding
            new_info=info.encode('UTF-8')
            d.write(new_info)
        except UnicodeEncodeError:
            num=len(info)
            #base condition
            if num==1:
                try:
                    d.write(info)
                except UnicodeEncodeError:
                    try:
                        new_info=info.encode('UTF-8')
                        d.write(new_info)
                    except UnicodeEncodeError:
                        d.write('*')
            else:
                writting(info[:num/2],d)
                writting(info[num/2:],d)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#write info into filename by mode 'w'
def writefile(info,filename):
    d=open(filename,'w')
    writting(info,d)
    d.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#write info into filename by mode 'w'
def append_file(info,filename):
    d=open(filename,'a')
    writting(info,d)
    d.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#just change the function name
#write info into filename by mode 'w'
def write_file(info,filename):
    d=open(filename,'w')
    writting(info,d)
    d.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#change list to dict, list looks like 'name  age sex email...'
#the first column is the index, like dict['name']=='age sex email...'
#nottab used to separate the content. if no, use '\t'
def list2dict(alist,*nottab):
    adict={}
    for ele in alist:
        item=ele.split()
        index=item[0]
        if len(nottab)==0:#no symbol for separation
            content='\t'.join(item[1:])
        else:
            sym=nottab[0]
            content=sym.join(item[1:])
        adict[index]=content
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if the file doesn't exit or size=0, wait until size>0
#after size, wait 1s to check whether the size of file change
def wait_file(fn,quiet=False,checksize_interval=1):
    #if file is not there
    if not filexist(fn):
        if not quiet:
            print('waiting the system to write the file named \''+fn+'\'...\n')
        a=time.time()
        count=1
        while (not filexist(fn)):
            b=time.time()
            #every 5 min, print once
            if b-a> 300:
                a=b
                minu=str(count*5)
                if not quiet:
                    print(minu+' min passed, still waiting...\n')
                count+=1
        ifchange='yes'
        while (ifchange=='yes'):
            statinfo=os.stat(fn)
            size1=statinfo.st_size
            time.sleep(checksize_interval)
            statinfo=os.stat(fn)
            size2=statinfo.st_size
            if size1==size2:
                ifchange='No'
    else:#if it is there, check file size
        ifchange='yes'
        while (ifchange=='yes'):
            statinfo=os.stat(fn)
            size1=statinfo.st_size
            time.sleep(checksize_interval)
            statinfo=os.stat(fn)
            size2=statinfo.st_size
            if size1==size2:
                ifchange='No'

    if not quiet:
         print('writing to '+fn+' finished! begin following process!\n')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#read the output of the cmd. the cmd should print output on the screen
def read_cmd(cmd):
    info=''
    f=os.popen(cmd)
    for line in f:
        info+=line
    f.close()
    return info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#read the output of the cmd. the cmd should print output on the screen
def cmd2info(cmd):
    info=''
    f=os.popen(cmd)
    for line in f:
        info+=line
    f.close()
    return info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#split the cmd output by '\n' into a list
def cmd2lines(cmd):
    astring=cmd2info(cmd).strip()
    if astring=='':#''.split('\n')will=['']
        return []
    else:
        return astring.split('\n')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#change dict to list. the dict looks like dict[index]=content
#the output list looks like: ['index\tcontent']
#nottab used to separate the content. if no, use '\t'
def dict2list(adict,*nottab):
    alist=[]
    for index in adict:
        if len(nottab)==0:#no symbol indicated in the input
            line=index+'\t'+adict[index]
        else:
            sym=nottab[0]
            line=index+sym+adict[index]
        alist+=[line]
    return alist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#convert dict[key]=[a,b] to 'key\ta,b\n'
#or adict[key]=a to 'key\ta\n'
def dict2info(adict):
    info=''
    for key in adict:
        value=adict[key]
        if type(value) is list:
            if len(value)==0:
                value='[]'
            elif type(value[0]) is list:
                value=str(value)
            else:#[1,2,3]
                value=','.join(value)
        if type(value) is dict:
            alist=[]
            for kkk in value:
                alist+=[str(kkk)+':'+str(value[kkk])]
            value='\t'.join(alist)

        info+="%s\t%s\n" % (key,value)
    return info

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#backward of info2dict
#the first as key, others as value
#you can indicate your syntax for spliting text
#if more than 2 items in one line, then the value is a list
#all list True if you want all values to be list
def info2dict(info,syntax=[],all_list=False):
    if type(syntax) is str:
        syntax=[syntax]
    lines=info.split('\n')
    lines=del_empty(lines)

    adict={}
    for line in lines:
        for syn in syntax:
            line=line.replace(syn,' ')
        item=line.split()
        key=item[0]
        if len(item)>=2:
            if len(item)==2:
                theitem=item[1]
                if theitem.startswith('[') and theitem.endswith(']'): #is a list
                    cmd='value=%s' % theitem
                    exec(cmd)
                elif all_list:
                    value=[theitem]
                else:
                    value=theitem
            else:
                value=item[1:]
            adict[key]=value
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def dict2line(adict):
    alist=[]
    for a in adict:
        alist+=['%s:%s' % (a,adict[a])]
    info='\t'.join(alist)
    return info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#submit form to the url (which directly accept the form), return respond
#form is a dict
def web_submit(form,url):
    import urllib.request, urllib.parse, urllib.error,urllib.request,urllib.error,urllib.parse
    data = urllib.parse.urlencode(form)
    request = urllib.request.Request(url, data)
    response = urllib.request.urlopen(request)
    page = response.read()
    return page
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#sort blist by the value in alist
#alist and blist has the, same length
#the order is from small to big
def sort_key(alist,blist,reverse=False):
    adict={}
    for i,key in enumerate(alist):
        if key not in adict:
            adict[key]=[blist[i]]
        else:
            adict[key]+=[blist[i]]
    alist=unify(alist)
    alist.sort()
    #print dict2info(adict)
    if reverse:
        alist.reverse()
    new_list=[]
    #count=0
    for key in alist:
        new_list+=adict[key]
    return new_list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#skip_inner_empty
def file2lines(fn,skip_inner_empty=False,skip_last_empty=True):
    info=txt_read(fn)
    if skip_inner_empty:#skip all empty lines
        lines=info.strip().split('\n')
        new_lines=[]
        for line in lines:
            if line.strip()!='':
                new_lines+=[line]
        return new_lines
    elif skip_last_empty:#only skip the last empty lines
        return info.strip().split('\n')
    else:
        return info.split('\n')

#change info into a stringIO object
def info2fp(info):
    import io
    return io.StringIO(info)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#copy file to my server directory
def gofigure(fn):
    os.system('cp %s /data/www/wenlin/html/figures/' % fn )

#use os.system to make dir
#if exist, do not make
def mkdir(wdir,new=False):
    if type(wdir) is list:
        for eachdir in wdir:
            mkdir(eachdir)
    elif type(wdir) is str:
        if not os.path.exists(wdir):
            os.system('mkdir -p %s' % wdir)
        elif not os.path.isdir(wdir):
            print('%s is not a directory!' %  wdir, file=sys.stderr)
        elif new==True:#if have the dir and need to delete the contents
            os.system('rm -r %s/*' % wdir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#overlap of two lists,regardless of element order
def overlap(alist,blist):
    overlap_list=[]
    if len(alist)>len(blist):
        biglist=alist
    else:
        biglist=blist
    #print alist,'\n',blist,'\n',clist
    #print biglist
    for item in biglist:
        if (item in alist) and (item in blist):
            overlap_list+=[item]
            #print item
    return overlap_list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#construct a dict like {'a':['b'],.....}
#if key not in dict, create one
def put_in_dict(adict,key,value):
    try:
        assert(adict[key])#if key in dict
        adict[key].append(value)
    except KeyError:
        adict[key]=[value]
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#check how many jobs in the queue, label is to identify jobs
def check_queue(label=None,rlabel=None):
    base='qstat -u "*" '
    if label==None:
        cmd=base
    elif type(label) is str:
        cmd='%s| grep "%s" ' % (base,label)
    elif type(label) is list:
        cmd=base
        for a in label:
            cmd='%s| grep "%s" ' % (cmd,a)

    if rlabel!=None:
        if type(rlabel) is str:
            cmd='%s| grep -v "%s" ' % (cmd,label)
        if type(rlabel) is list:
            for a in rlabel:
                cmd='%s| grep -v "%s" ' % (cmd,a)

    number=int(cmd2info('%s |wc -l' % cmd))
    return number
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#used to control the job in queue
#label is to identify jobs, number is the maxiumal allow, wait is the waiting time (second)
#if job more than number, wait
def watch_job(label="wenlin",rlabel=None,number=30,wait=200):
    job_number=check_queue(label,rlabel)
    while(job_number>=number):#while overload, wait
        time.sleep(wait)
        job_number=check_queue(label,rlabel)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#replace the chars between start_label and end_label
#if b == None, take until the end of the line
#if a == None, take from the start
def replace_between(string, repl, start_label=None, end_label=None):
    i = 0
    j = len(string)
    if start_label:
        i = string.find(start_label) + len(start_label)
    if end_label:
        j = string.find(end_label, i)

    if i == -1:
        print("label %s is not found!" % start_label, file=sys.stderr)
        return ""
    if j == -1:
        print("label %s is not found!" % end_label, file=sys.stderr)
        return ""
    newstring = "%s%s%s" % (string[:i], repl, string[j:])
    return newstring
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#take the chars between start_label and end_label
#if b == None, take until the end of the line
#if a == None, take from the start
def find_between(string, start_label=None, end_label=None):
    i = 0
    j = len(string)
    if start_label:
        i = string.find(start_label) + len(start_label)
    if end_label:
        j = string.find(end_label, i)

    if i == -1:
        print("label %s is not found in %s!" % (start_label, string), file=sys.stderr)
        return ""
    if j == -1:
        print("label %s is not found in %s!" % (start_label, string), file=sys.stderr)
        return ""
    return string[i:j]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pickle_write(obj, fn):
    f = open(fn, "wb")
    pickle.dump(obj, f)
    f.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pickle_read(fn):
    f = open(fn, "rb")
    obj = pickle.load(f)
    f.close()
    return obj
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Caution: only applicable to adict[a] = b.(1 to 1 mapping)
def reverse_dict(adict):
    bdict = {}
    for key in adict:
        bdict[adict[key]] = key
    return bdict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ordered_dict(object):
    """
    a class with value in order
    """
    def __init__(self, adict):
        self.adict = adict
        for key in adict:
            if type(key) is int:
                raise TypeError("int type key for dictionary is not accepted")
        self.order_list = self._make_order(adict)

    def _make_order(self, adict):
        """order the dict by the value"""
        keys = []
        values = []
        for key in adict:
            keys.append(key)
            values.append(adict[key])
        keys = sort_key(values, keys)
        final_list = []
        for key in keys:
            final_list.append([key, adict[key]])
        return final_list

    def __getitem__(self, index):
        return self.order_list[index]

def csv_read(fn):
    f = open(fn, 'rb')
    freader = csv.reader(f)
    alist = []
    for row in freader:
        alist.append(row)
    f.close()
    return alist

def csv_write(listoflist, dn):
    f = open(dn, 'wb')
    fwriter = csv.writer(f)
    for row in listoflist:
        fwriter.writerow(row)
    f.close()

def run(cmd):
    return os.system(cmd)

def write_lines(lines, dn):
    write_file('\n'.join(lines), dn)


def lastName(fn):
    fn = fn.strip('/ ').split('/')[-1]
    return fn
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def voidi():
    return 1








