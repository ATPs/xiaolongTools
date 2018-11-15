#lib to deal with lab frame using bootstrip
#docs see: http://twitter.github.com/bootstrap/getting-started.html

#used in formating table
from .basic import basic_table as bt
from .basic import wrap_page

class general(object):
    """
    general class for all bootstrap classes.
    includes the lib location.
    """
    def __init__(self):
        self.inhead=self.make_inhead()

    def make_inhead(self):
        """
        the libs required for running bootstrap
        """

        a = '''
        <link href="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/css/bootstrap.css" rel="stylesheet" media="screen">
        <link href="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/css/docs.css" rel="stylesheet" media="screen">
        <link href="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/css/bootstrap-responsive.css" rel="stylesheet" media="screen">
        <script src="http://code.jquery.com/jquery.js"></script>
        <script src="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/js/bootstrap.js"></script>
        '''

        #jquery location
        #a = '<script type="text/javascript" src="http://prodata.swmed.edu/wenlin/necessary_tools/jquery-1.8.2.min.js"></script>\n'
        #css table
        #a+= '<link rel=stylesheet type=text/css href="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/css/bootstrap.css">\n'
        #javascript lib
        #a+= '<script type="text/javascript" src="http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/js/bootstrap.js"></script>\n'
        #some image... I don't know the funciton yet
        #a+='http://prodata.swmed.edu/wenlin/necessary_tools/bootstrap/img/glyphicons-halflings.png'
        return a

    def mytable(self,alist,more='class="table table-striped table-bordered"'):
        """
        created bordered striped table with header
        """
        return bt(alist,more,with_head=True)

    def checkout(self,html='',title='Wenlin'):
        return wrap_page(html,inhead=self.inhead,title=title)


class homepage(general):
    """
    the class to building an easy homepage
    """
    def __init__(self):
        general.__init__(self)
        self.html=self.make_html(self,title="^_^")

    def make_html(self):
        pass
