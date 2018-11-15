#advanced action for shelve
#wrapped actions

import os,sys
shelve_path="/usr6/local/lib/python2.7/"
if shelve_path not in sys.path:
    sys.path.insert(0,shelve_path)


import shelve

class shelveIO(object):
    """
    wrapped read and write actions as if the object is a dict
    in the get_value mode, if no key in the dict, return None
    in the update mode, overwrite the old key and value
    """
    def __init__(self,fdict):
        #the file name of the dict
        self.fdict=fdict

    def get_value(self,key):
        d=shelve.open(self.fdict)
        try:
            value=d[key]
        except KeyError:
            value=None
        d.close()
        return value

    def update(self,newdict):
        d=shelve.open(self.fdict, writeback=True)
        for key in newdict:
            d[key]=newdict[key]
        d.close()
        os.system('chmod a+w %s' % self.fdict)

    def keys(self):
        d=shelve.open(self.fdict)
        a=list(d.keys())
        d.close()
        return a

