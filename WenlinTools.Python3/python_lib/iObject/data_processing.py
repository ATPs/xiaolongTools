#input a list of number,do something
#put the data into bin
#test only on int dataset, some problem exist for float type
class make_bin(object):
    def __init__(self,alist,spacing=1,start=None,end=None):
        self.keys,self.values,self.adict=self.make_result(alist,start,end,spacing)

    def make_result(self,alist,start,end,spacing):
        #make the bin
        if start==None:
            start=int(min(alist))
        if end == None:
            end=int(max(alist))+spacing

        keylist=[]
        keylist.append('<%s' % start)

        adict={}
        adict['<%s' % start]=0
        adict['>%s' % end]=0

        i=start
        while(i<end):
            adict['%s<x<%s' % (i,i+spacing)]=0
            keylist.append('%s<x<%s' % (i,i+spacing))
            i += spacing

        keylist.append('>%s' % end)


        for element in alist:
            placing=int(round((element-start)/spacing,0))
            if element<start:
                adict['<%s' % start]+=1
            elif element>end:
                adict['>%s' % end]+=1
            else:
                i=placing*spacing+start
                if i==end:#if round make it bigger
                    adict['%s<x<%s' % (end-spacing,end)]+=1
                else:
                    try:
                        adict['%s<x<%s' % (i,i+spacing)]+=1
                    except:
                        print('error for element %s' % element)

        values=[]
        for key in keylist:
            values.append(adict[key])
        return keylist,values,adict

    #user-friendly print out
    def checkout(self):
        alist=[]
        for key in self.keys:
            alist.append('%s\t%s'%(key,self.adict[key]))
        return '\n'.join(alist)




