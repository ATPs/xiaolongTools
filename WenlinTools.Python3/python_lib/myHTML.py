#class to build the HTML webpage


#input:
#output:
#algorithm:
#author:wenlin; Date:2012-


#functions

#add <div> tag
#def add_div(info, more=''):
#    info_n='<div '+more+' >\n'+info+'</div>\n'
#    return info_n

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class basicHTML:
    def __init__(self,title='Wenlin',tag='basic'):
        self.start=self.create_start(title,tag)
        self.end=self.create_end()
        self.content=''


    def create_start(self,title,tag):
        info_start='''<!-- template borrowed from Qian -->
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<meta name="keywords" content="content_test"/>
<title>title_test</title>
<link rel=stylesheet type=text/css href=http://prodata.swmed.edu/wenlin/necessary_tools/styleSEQ.css>
<script src="http://prodata.swmed.edu/wenlin/necessary_tools/jmol/Jmol.js" type="text/javascript"></script>
<script src="http://prodata.swmed.edu/wenlin/necessary_tools/javascripts/ajaxgold.js" type="text/javascript"></script>
<script src="http://prodata.swmed.edu/wenlin/necessary_tools/javascripts/utils.js" type="text/javascript"></script>
<script src="http://prodata.swmed.edu/wenlin/necessary_tools/javascripts/myUtils.js" type="text/javascript"></script>
</head>
<body>
'''
        info=info_start.replace('title_test',title)
        if tag=='server':
            info=info.replace('styleSEQ.css','server.css')    
        return info
    
    def body_bg(self,more=''):
        if more=='':
            #default
            bg='<body background=info/promals_bg.gif style="background-color: #E6F2FF; position: relative; text-align: center;">'
        else:
            bg=more
        self.start=self.start.replace('<body>',bg)

    def create_end(self):
        info_end='''
</body>
</html>'''
        return info_end

    def header(self, head='Wenlin',t='h2', more=''):
        info_head='''<div '''+more+'''>
<h1>header</h1>
</div>
'''
        info_head=info_head.replace('header',head)
        info_head=info_head.replace('h1',t)
        self.content+=info_head

    def add_div(self,info, more=''):
        if more!='':
            info='<div '+more+' >\n'+info+'</div>\n'
        else:
            info='<div>\n'+info+'</div>\n'
        return info

    aaaaa='''
    #wrap each section for the server style
    def section_wrap(self, title, info, more='class="divisions"'):
        more=' '+more+' '
        title='<div class="divisionTitle">\n'+title+'\n</div>\n'
        info=title+info
        info='<div'+more+'>\n'+info
        info+='</div>\n'
        return info
'''
    def section_wrap(self,title='',content=''):
        aa='''<div class="divisionTitle">
        empty_section 
        </div>
        <p class="param">
        moreContent
        </p>
        '''
        if title!='':
            aa=aa.replace('empty_section',title)
        if content!='':
            aa=aa.replace('moreContent',content)
        aa=self.add_div(aa, 'class="divisions"')
        return aa



    #here info=alist, because I copy these from table2html.py
    def easytable(self, info):
        count=0#0=has count, 1=no count

        table1='<div>\n<table class="easy" border="1">\n'

        for line in info:
            if line!='':
                table1+='<tr>\n'
                line=line.split('\t')
                for element in line:
                    if count==0:#count line
                        table1+='<th class="headcell">'+element+'</th>\n'
                    elif count%2==0:#the even line
                        table1+='<td class="cell1">'+element+'</td>\n'
                    else:
                        table1+='<td>'+element+'</td>\n'
                table1+='</tr>\n'
                count+=1
        table1+='</table>\n</div>\n'
        self.content+=table1

    def linktable(self, alist):
        pass

    #use 'addform\n' as a tag to add content
    #when finish form, use form_finish to remove tag
    def form_frame(self, action, style=''):
        frame='<form enctype="multipart/form-data" action="'+action+'" method="POST"'
        if style=='':
            frame+=' class="basic_form" >\naddform\n</form>\n'
        else:
            frame+=' class="'+style+'">\naddform\n</form>\n'
        frame=self.add_div(frame)
        self.content+=frame
    
    def add_form(self,info):
        a=self.content
        self.content=a.replace('addform',info+'addform')
    
    def form_finish(self):
        a=self.content
        self.content=a.replace('addform\n','')

    def form_submit(self,unknown):
        pass

    def hyperlink(self, word, link, in_sentence=''):
        linkformat='<a href="javascript:void(0)" onclick=\'openwindow("'+link+'")\'>'+word+'</a>'
        if in_sentence=='':
            self.content=self.content.replace(word, linkformat)
        else:#indicate the sentence
            linksentence=in_sentence.replace(word, linkformat)
            self.content=self.content.replace(in_sentence,linksentence)

    def write(self, filename):
        info=self.start+self.content+self.end
        if '.html' not in filename:
            filename+='.html'
        d=open(filename,'w')
        d.write(info)
        d.close()
    
    def page(self):
        info=self.start+self.content+self.end
        return info 



#used to generate well-format template
class toggle(basicHTML):
    def __init__(self):
        pass

    def upload_fasta(self, title='text input',more=''):
        aa='''<p style="position: relative; margin-top: 0px;margin-down:2px;" class="optitle">
        <b>%s:</b>
        <input type="button" style="position: absolute; right: 12px; top: -6px; font-size: 11px;" onclick="clear_query()" value="Clear text input">''' % title

        aa+='''</p>
        <div style="width: 100%;"><br>
        <textarea style="width: 100%;" rows="10" name="upload_text"></textarea><br>
        </div><br>
        Or upload a file <input type="file" size="24" name="upload_file">
        <br>
        ''' 
        aa=self.add_div(aa,more)
        self.upload_fasta_data=aa
        return aa

    def submit(self, more=''):
        aa='''<div class="divisionTitle">
        DATA SUBMIT 
        </div>
        <p class="param">
        Enter an <b>email</b> to receive the result (required): <input type="text" size="20" name="email"><br><br>
        Enter a <b>job name</b> (recommended):
        <input type="text" size="10" name="job_name">
        <input type="submit" class="btn" value="Submit"><input type="reset" class="btc" value="Reset">

        </p>
        '''
        aa=self.add_div(aa,'class="divisions"')
        self.submit_data=aa
        return aa

    def para_wrap(self,content,title='PARAMETERS'):
        aa='<div class="divisionTitle" style="width: 200px">\n%s\n' % title
        aa+='<a onclick=\'toggle_region("params_all")\' href="javascript:void(0)" style="text-decoration">show/hide</a>\n'
        aa+='</div>\n'
        
        aa+='<div id="params_all" style="display: none;">\n'
        aa+=content
        aa+='</div>\n'
        return aa
        

    def para_list(self,adict,more=''): 
        if more=='':
            bb='<ul class="param">\n'
        else:
            bb='<ul %s>\n' % more
        
        for key in adict:
            bb+='<li>\n'
            bb+=key+':'
            cc=adict[key]
            items=cc.split('=')
            bb+='<input type="text" size="4" value="%s" name="para_%s">\n' % (items[1],items[0])
            bb+='</li>\n'
        bb+='</ul>'

        return bb

    def para_checkbox(self,adict,more=''):
        pass

    def para_radio(self,adict,more=''):
        pass

    #below is used in showing alignment
    
    ######
    def getmax(self,alist):
        win=0
        for adict in alist:
            a=int(adict['length'])
            if a>win:
                win=a
        return win        

    #alist like [adict, bdict], adict like {qb:20,qe:70,begin:1,end:50,length:90}
    #the first is query. adjusted range according to query
    #return a list of prepared element
    def cover_bar(self,alist, width='400px'):
        window=int(width.replace('px',''))
        maxlen=self.getmax(alist)
        querylen=int(alist[0]['length'])
        if querylen>0.6*window and maxlen>window: 
            ratio=querylen/(0.6*window)
        else:
            ratio=1

        print('ratio is %s\n' % ratio)    
        qr=querylen/ratio
        #every span consist of 1,filler(white),2,noalileft([==),3,align;4,noaliright    
        bar_list=[]
        for i,aln in enumerate(alist):
            if i==0:#the query
                qstart=(int(aln['begin'])-1)/ratio
                qend=int(aln['end'])/ratio
                qlen=int(aln['length'])/ratio
                each_end=(window-qlen)/2
                info='lfiller:%s\talign:%s\trfiller:%s' % (qstart+each_end,qend-qstart,window-qend-each_end)
                tt1=self.draw_bar(info,window)
                tt='<nobr>\n%s</nobr>\n' % tt1
                bar_list+=[tt]
            else:
                start=(int(aln['begin'])-1)/ratio
                end=int(aln['end'])/ratio
                length=int(aln['length'])/ratio
                qstart=int(aln['qb'])/ratio
                qend=int(aln['qe'])/ratio
                shift=qstart-start+each_end
                a=[shift,start,end-start,length-end,window-length-shift]
                info='lfiller:%s\tnoalileft:%s\talign:%s\tnoaliright:%s\trfiller:%s' % (a[0],a[1],a[2],a[3],a[4])
                tt1=self.draw_bar(info,window)
                tt='<nobr>\n%s</nobr>\n' % tt1
                bar_list+=[tt]
        return bar_list

    def draw_bar(self,info,window):
        #transfer to a dict
        items=info.split('\t')
        adict={}
        for item in items:
            vv=item.split(':')
            adict[vv[0]]=vv[1]
        print(adict)    
        #now really drawing in the order of:
        #filler,noalileft,align,noalight,filler
        bar_info=''
        if 'noalileft' not in adict:#if it is the query
            bar_info+=self.draw_span('filler',adict['lfiller'])+self.draw_span('align',adict['align'])+self.draw_span('filler',adict['rfiller'])+'\n'
        else:
            #add into bar_info one by one
            #1,deal with left part
            if int(adict['lfiller'])>0:#if <0, do not add it
                bar_info+=self.draw_span('filler',adict['lfiller'])+self.draw_span('noalileft',adict['noalileft'])
            else:#now lfiller is the shifted value
                shift=adict['lfiller']
                bar_info+=self.draw_span('noalileft','40')+self.draw_span('filler','20',rs=str(0-int(shift)-40))+self.draw_span('noaliright',str(int(adict['noalileft'])+int(shift)-60))
            
            #2,middle
            bar_info+=self.draw_span('align',adict['align'])
            
            #3,right part
            if int(adict['rfiller'])>0:#if length not exceed
                bar_info+=self.draw_span('noaliright',adict['noaliright'])+self.draw_span('filler',adict['rfiller'])
            else:#now lfiller is the shifted value
                shift=adict['rfiller']
                bar_info+=self.draw_span('noalileft',str(int(adict['noaliright'])+int(shift)-60))+self.draw_span('filler','20',rs=str(0-40-int(shift)))+self.draw_span('noaliright','40')
           
        return bar_info

    def draw_span(self, key, lg, rs=''):
        if rs=='':
            info='<span class="%s" style="width:%spx"></span>' % (key,lg)   
        else:
            info='<span class="%s" title="%s residues" style="width:%spx"></span>\n' % (key,rs,lg)
        return info    



#define diffent component as toggle
#combine the element hierarchily by different type        
#use form first as example to develop the toggle structure(2012-2-4)
class advancedHTML(basicHTML):

    #combine different toggle by type
    #this is the basic function to build the hierarchy
    def combine(self, a , b, ltype=''):    
        if a[-1]!='\n':
            content=a+'\n'+b
        else:
            content=a+b
        
        #no label indicated
        if ltype=='':
            content=self.add_div()
            return content











    









