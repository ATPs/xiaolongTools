#general function for all the chart type
import numpy

class general(object):
    def __init__(self,title="^_^"):
        self.title=title
        self.inhead=self.make_inhead(title)
        self.javascript=''

    def make_inhead(self,title):
        a = """
    <title>%s</title>

    <meta http-equiv="X-UA-Compatible" content="chrome=1">
    <meta http-equiv="CACHE-CONTROL" CONTENT="NO-CACHE">
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" >
    <meta name="keywords" content="canvasxpress, canvas, html5, graph, chart, plot, javascript, javascript library, genomic, scientific, android, animation, bar graph, line graph, dotplot, boxplot, heatmap, newick, scatter, 3d, pie, correlation, venn, network, market, candlestick, genome browser, isaac neuhaus">
    <meta name="description" content="">
    <meta http-equiv="Content-Language" content="en-us" >

    <script type="text/javascript" src="http://prodata.swmed.edu/wenlin/necessary_tools/canvasXpress-SF/js/canvasXpress.min.js"></script>
	""" % title
        return a

    def checkout(self):
        a = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
    <html>
    """
        a += '<head>\n%s<script>var showing = function (){\n%s\n}</script>\n</head>\n' % (self.inhead, self.javascript)
        a += '<body><div style="margin: auto;width: 928px;"><h2 align="center">%s</h2><canvas id="canvas" width="613" height="500"></canvas></div><script>showing()</script></body>\n</html>\n' % self.title
        return a

class basic_heatmap(general):
    """
    a class building the framework of other heatmap classes
    x,z is for x and y axis
    for more info, refer to http://canvasxpress.org/documentation.html
    """
    def __init__(self,title='^_^',data={},x='{}',z='{}',t='{}',config={},events='{}'):
        general.__init__(self,title)
        #data in json format, equal to "y" in online tutorial
        self.y='y: %s' % data
        #the label of y axis
        self.x='x: %s' % x
        #the label of x axis
        self.z='z: %s' % z
        #the distance tree
        self.t='t: %s' % t

        self.config=config
        self.events=events

    def checkout(self,color_grad="white-blue"):
        config={"graphType": "Heatmap","heatmapType": color_grad,"varLabelRotate": 45,"overlaysWidth": 20, }
        config.update(self.config)
        self.javascript="new CanvasXpress('canvas', {%s, %s, %s, %s}, %s, %s )" % (self.x, self.y, self.z, self.t, str(config), self.events)
        return super(basic_heatmap,self).checkout()


class interact_map(basic_heatmap):
    #the index in interaction begins with 0
    #the interaction is {'i j':number,}
    def __init__(self,title="^_^", length=0, interaction={},threshold=50.7758521,label=None):
        if label!=None:
            self.label=label
        else:
            self.label="distance"
        data=self.make_data(length,interaction,threshold)
        basic_heatmap.__init__(self,title,data)

    def make_data(self,length,interaction,threshold):
        data={}
        #print interaction
        axis=[str(i+1) for i in range(length)]
        data["vars"]=axis
        data["smps"]=axis
        data["desc"]=[self.label]
        data["data"]=[]
        for i in range(length):
            a1=[]
            for j in range(length):
                key1='%s %s' % (i,j)
                try:
                    value=interaction[key1]
                    #inforce the threshold
                    if value>threshold:
                        value=threshold
                except:
                    value=threshold
                a1+=[value]
            data["data"]+=[a1]
        return data

class easy_heatmap(basic_heatmap):
    """
    input: x,y value and a matrix
    """
    def __init__(self,title="^_^",xAxis=[],yAxis=[],z={},label=None,threshold=None):
        if label!=None:
            self.label=label
        else:
            self.label="Grad"


        data=self.make_data(xAxis,yAxis,z,threshold)
        basic_heatmap.__init__(self,title,data)


    def make_data(self,xAxis=[],yAxis=[],z={},threshold=None):
        data={}
        data["vars"]=xAxis
        data["smps"]=yAxis
        data["desc"]=[self.label]
        data["data"]=[]
        #need to test if M!=N for matrix
        if isinstance(z,numpy.ndarray):#if it is a matrix
            M,N=z.shape
            if data["vars"]==[]:
                data["vars"]=[i+1 for i in range(M)]
            if data["smps"]==[]:
                data["smps"]=[i+1 for i in range(N)]
            print(z)
            for i in range(M):
                a1=[]
                for j in range(N):
                    point=z[i,j]
                    if threshold!=None and point<threshold:
                        a1+=[z[i,j]]
                    elif threshold!=None:
                        a1+=[threshold]
                    else:
                        a1+=[point]
                data["data"]+=[a1]
        else:
            print('not numpy.ndarray instance, need more work here!')
        return data

