#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#start paralleled threads, wait till all finished
#function_list=[func1,func2]
#argv_list=[[arg1,arg2],[arg1,arg2]]
def old_run_threads(function_list,argv_list):
    from threading import Thread
    i=-1
    pointers=[]
    for function in function_list:
        i+=1
        cmd='t%s=Thread(target=function,args=%s)' % (i,argv_list[i])
        exec(cmd)
        exec('t%s.start()' % i)
        pointers+=['t%s' % i]
    #wait all finished
    for p in pointers:
        exec('%s.join()' % p)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#from '(' to the end
def get_argv(line):
    argv=''
    start=False
    for char in line:
        if char=='(':
            start=True
        if start:
            argv+=char
    argv='['+argv[1:-1]+']'
    name=line.split('(')[0]
    return name,argv            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#start paralleled threads, wait till all finished
#function_list=['func1(arg1,arg2)','func2(arg1,arg2)']
def run_threads(function_list):
    from threading import Thread
    pointers=[]
    i=-1
    for function in function_list:
        i+=1
        name,argv=get_argv(function)
        cmd='t%s=Thread(target=%s,args=%s)' % (i,name,argv)
        print(cmd)
        exec(cmd)
        exec('t%s.start()' % i)
        pointers+=['t%s' % i]
    #wait all finished
    for p in pointers:
        exec('%s.join()' % p)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





