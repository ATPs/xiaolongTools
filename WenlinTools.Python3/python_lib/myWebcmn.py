import cgi,cgitb,os,sys
cgitb.enable()

#two environment valuable to run qsub
os.environ["SGE_ROOT"] = "/opt/gridengine"
os.environ["SGE_QMASTER_PORT"] = "536"
os.environ["SGE_EXECD_PORT"] = "537"

def print_config():
    #standard format for html
    print("Content-Type: text/html")
    print('')

def path_config(sys):
    paths=[]
    #python path
    paths.append('/home/wenlin/my_programs/python_lib')
    #biopython path
    paths.append('/home/wenlin/local/biopython-1.58')
    sys.path+=paths
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

def basic_table(alist, more=""):
    """
    basic function to parse a web table out. all the style options can be changed.
    the input structure for alist is like [['a','b'],['c','d']], 'more' is for the overall table style.
    if you want to specify the style for some row, then change that row to be a list with two element, like [['style="center"', ['a','b'] ] ,['c','d']]
    if you want to specify the style for some column, then change that column to be a list, like [[['style="center"','a'],'b'],['c','d']]
    """
    info='<table %s>\n' % more
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
                info+='<td %s>%s</td>\n' % tuple(column)
            else:
                info+='<td>%s</td>\n' % column
        info+='</tr>\n'
    info+='</table>\n'
    return info


def wrap_page(info,title='^_^',inhead=None, infoot=None):
    a='<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">\n'
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


def sendMail(subject="Wenlin",content="Hello",sender="wenlin.li@utsw.edu",receiver=""):
    import smtplib

    if '@' not in sender:
        fromaddr = "%s@prodata.swmed.edu" % sender
    else:
        fromaddr = sender
    toaddr = receiver

    msg = ("From: %s\r\nTo: %s\r\nSubject: %s\r\n" % (fromaddr, toaddr, subject))
    msg += content
    server = smtplib.SMTP('smtp.swmed.edu','25')
    server.sendmail(fromaddr, toaddr, msg)
    server.quit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_bootstrap():
    info = '<link type="text/css" rel="stylesheet" href="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/css/bootstrap.min.css">'
    info += '<script type="text/javascript" src="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/js/bootstrap.min.js"></script>'
    return info

def get_bootstrap2():
    info = '<link type="text/css" rel="stylesheet" href="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap-3.1.0/css/bootstrap.min.css">'
    info += '<script type="text/javascript" src="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap-3.1.0/js/bootstrap.min.js"></script>'
    return info

def hyperlink(text, link, title=None, new=True):
    if title == None:
        html = '<a target="_blank" href="%s">%s</a>' % (link, text)
    else:
        html = '<a target="_blank" href="%s" title="%s">%s</a>' % (link, title, text)
    if not new:
        html = html.replace('target="_blank" ', '')
    return html



def check_input(string, allow=[]):
    """
    only allow: 1. number 2. digit 3. dot 4. dashes 5. underscore
    if need more allowance, put in allow
    """
    allowed = set([".","_","-"]) | set(allow)
    flag = True
    for char in string:
        if not char.isalnum():
            if char not in allowed:
                flag = False
    return flag



def voi():
    pass

