#copy from Webcmn.py
#only contant the function and class

import cgi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def transform():
    form=cgi.FieldStorage()
    adict={}
    for key in form:
        adict[key]=form.getvalue(key,'')
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def clear_page(content=''):
    print("""<script type="text/javascript">
                document.body.innerHTML = '%s';
                        </script>""" % content)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#used in the next section
import urllib.request, urllib.parse, urllib.error,urllib.request,urllib.error,urllib.parse
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def url_request(url,form=None):
    #submit a url request
    #a form submission can be included

    #header={'User-Agent':'Mozilla/5.0 (Windows; U; Windows NT 6.1; en-US; rv:1.9.1.6) Gecko/20091201 Firefox/3.5.6'}
    if form!=None:
        data = urllib.parse.urlencode(form)
        request = urllib.request.Request(url, data)
    else:
        request = urllib.request.Request(url)
    response = urllib.request.urlopen(request)
    page = response.read()
    #type = sys.getfilesystemencoding()
    #page = page.decode("UTF-8").encode('gbk')
    return page

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def basic_table(alist, more="", with_head=False):
    """
    basic function to parse a web table out. all the style options can be changed.
    the input structure for alist is like [['a','b'],['c','d']], 'more' is for the overall table style.
    if you want to specify the style for some row, then change that row to be a list with two element, like [['style="center"', ['a','b'] ] ,['c','d']]
    if you want to specify the style for some column, then change that column to be a list, like [[['style="center"','a'],'b'],['c','d']]
    if header included, set with_head to be true
    """
    info='<table %s>\n' % more
    if with_head:#if header present in the list
        info+='<thead>\n'
        row=alist[0]
        alist=alist[1:]
        #copied from below

        #check whether ['a','b'] or ['style="sth."',['a','b']]
        if type(row[-1]) is list:
            info+='<tr %s>\n' % row[0]
            row=row[-1]
        else:
            info+='<tr>\n'
        for column in row:
            #check whether 'a' or ['style="sth."','a']
            if type(column) is list:
                info+='<th %s>%s</th>\n' % tuple(column)
            else:
                info+='<th>%s</th>\n' % column

        info+='</thead>\n'
        info+='<tbody>'

    for row in alist:
        #check whether ['a','b'] or ['style="sth."',['a','b']]
        if type(row[-1]) is list:
            info+='<tr %s>\n' % row[0]
            row=row[-1]
        else:
            info+='<tr>\n'
        for column in row:
            #check whether 'a' or ['style="sth."','a']
            if type(column) is list:
                info+='<th %s>%s</th>\n' % tuple(column)
            else:
                info+='<th>%s</th>\n' % column
        info+='</tr>\n'

    if with_head: info+='</tbody>'

    info+='</table>\n'
    return info
def wrap_page(info,title='^_^',inhead=None, infoot=None):
    a='<!DOCTYPE html>\n'
    a+='<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">\n'
    a+='<head>\n'
    a+='<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" >\n'
    a+='<title>%s</title>\n' % title
    a+='<link rel=stylesheet type=text/css href="http://prodata.swmed.edu/wenlin/necessary_tools/styleSEQ.css">\n'
    a+='<script src="http://prodata.swmed.edu/wenlin/necessary_tools/javascripts/ajaxgold.js" type="text/javascript"></script>\n'
    a+='<script src="http://prodata.swmed.edu/wenlin/necessary_tools/javascripts/utils.js" type="text/javascript"></script>\n'

    if inhead!=None:
        if type(inhead) is str:
            a+='%s\n' % inhead
        elif type(inhead) is list:# a list of lib in head
            for ele in inhead:
                a+='%s\n' % ele

    a+='</head>\n<body>\n'

    if infoot!=None:
        if type(infoot) is str:
            a+='%s\n' % inhead
        elif type(infoot) is list:# a list of lib in head
            for ele in infoot:
                a+='%s\n' % ele

    a+='%s\n</body>\n</html>\n' % info
    return a

def add_div(html,more=None):
    if more:
        return '<div %s>\n%s</div>\n' % (more,html)
    else:
        return '<div>\n%s</div>\n' % (html)



#a decorator to add <div> to text output
a="""
def add_div(style=None):
    def decorator(fn):
        def wrapper():
            if style:
                return '<div %s>\n%s\n</div>\n' % (style,fn())
            else:
                return '<div>\n%s</div>\n' % fn()
        return wrapper
    return decorator
"""
def voi():
    pass

