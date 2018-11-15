#!/usr/bin/env python

#function:map DCA result (index begin with 0) to pdb (index with 1)
#function:also calculate the distance of the tops (incooperate from ??)
#caution: the sequence of the pdb and dca should be the same
#input:dca result, pdb file, top number
#output: pml file with mapped interaction
#algorithm: reorder the index of pdb file, then mapping
#author:wenlin; Date:2012-8-24

#this file is modified based on evaluate_DCA_on_PDB

import math

#used by PDB.interact_map()
from iHTML.canvasXpress import interact_map as itmap

#used by pdb2constraint
from iObject.parser import horiz as hz

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def parse_dca_result(fn):
    lines=cmn.txt_read(fn).strip().split('\n')
    alist=[]
    for line in lines:
        items=line.split()
        alist+=[[str(int(items[2])+1),str(int(items[3])+1)]]
    return alist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def put_in_dict(index,adict,data):
    if type(data) is dict:
        if index not in adict:
            adict[index]=data
        else:
            adict[index].update(data)
    else:
        if index in adict:
            adict[index]+=[data]
        else:
            adict[index]=[data]
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def shortest_distance(alist,blist):
    short=99999
    for i in alist:
        for j in blist:
            dist=distance(i,j)
            if dist<short:
                short=dist
    return short
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def color_top(pairlist,top=10):
    color_list=['red','orange','yellow','limon','green','cyan','blue','magenta','pink','wheat']
    k=len(color_list)
    cmd=''
    top=min(len(pairlist),top)
    for i in range(top):
        cmd+="cmd.color('%s','resi %s+%s')\n" % (color_list[i%k],pairlist[i][0],pairlist[i][-1])
    return cmd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def select_top(pairlist,top=10):
    cmd=''
    top=min(len(pairlist),top)
    for i in range(top):
        cmd+="cmd.select('pair%s','resi %s+%s')\n" % ((i+1),pairlist[i][0],pairlist[i][-1])
        cmd+="cmd.disable('pair%s')\n" % (i+1)
    return cmd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calculate_pair_distance(pair_list,dist_info,top=30):
    #parse the distance to be indexed
    lines=dist_info.strip().split('\n')
    adict={}
    for line in lines:
        items=line.split()
        if int(items[-2])<int(items[-1]):
            key='%s,%s' % (int(items[-2]),int(items[-1]))
        else:
            key='%s,%s' % (int(items[-1]),int(items[-2]))
        #print key
        adict[key]=line

    newline=[]
    for pair in pair_list[:top]:
        key='%s,%s' % (pair[0],pair[1])
        newline+=[adict[key]]
    return '\n'.join(newline)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def residue_one2three_dict():
    adict={}
    adict.update({'A':'ALA'})
    adict.update({'C':'CYS'})
    adict.update({'D':'ASP'})
    adict.update({'E':'GLU'})
    adict.update({'F':'PHE'})
    adict.update({'G':'GLY'})
    adict.update({'H':'HIS'})
    adict.update({'I':'ILE'})
    adict.update({'K':'LYS'})
    adict.update({'L':'LEU'})
    adict.update({'M':'MET'})
    adict.update({'N':'ASN'})
    adict.update({'P':'PRO'})
    adict.update({'Q':'GLN'})
    adict.update({'R':'ARG'})
    adict.update({'S':'SER'})
    adict.update({'T':'THR'})
    adict.update({'V':'VAL'})
    adict.update({'W':'TRP'})
    adict.update({'Y':'TYR'})
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################################################
#############all the function above is the##################
#############one remained in the old file ##################
############################################################

def residue_three2one_dict():
    adict=residue_one2three_dict()
    newdict={}
    for key in adict:
        newdict[adict[key]]=key
    newdict['ABA']='A'
    return newdict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os,sys
import cmn

#used by pdb2constraints
import numpy as np
#from csb.numeric import dihedral_angle as dan

class PDB(object):
    def __init__(self,info,quiet=True):
        if os.path.exists(info):#check whether it is a file name
            self.info=cmn.txt_read(info)
            input_type='file'
        else:
            self.info=info
            input_type='info'

        if not quiet:
            print('input file type is %s' % input_type, file=sys.stderr)

    def reindex(self,template=None):
        #make the index beginning from 1 and continuous
        #if template indicated, use the plain text sequence to adjust index
        lines=self.info.split('\n')
        newline=[]
        count=0
        i=0
        k=len(lines)
        line='yes'
        if template!=None:
            rdict=residue_three2one_dict()
            template=template.strip()
            kk=len(template)
        appeared_num=-99
        while (i<k):
            line=lines[i]
            if line.startswith('ATOM') or line.startswith('HETATM'):#check naming
                num=int(line[22:26].strip())
                if num!=appeared_num:#new atom start
                    appeared_num=num
                    count+=1
                    if template!=None and count<=kk:
                        pdb_resi=rdict[line[17:20].strip()]
                        seq_resi=template[count-1]
                        while(pdb_resi!=seq_resi):
                            count+=1
                            seq_resi=template[count-1]
                        print(count)
                #fill=' '*(4-len(str(count)))+str(count)#reindex
                fill='{:>4}'.format(count)
                newline+=['%s%s%s' % (line[:22],fill,line[26:])]
            elif line.startswith('TER'):#the chain ends
                #fill=' '*(4-len(str(count)))+str(count)#reindex
                fill='{:>4}'.format(count)
                newline+=['%s%s%s' % (line[:22],fill,line[26:])]
                if int(lines[i+1][22:26])==1:
                    count=0
            else:
                newline+=[line]
            i+=1
        self.info='\n'.join(newline)

    #make the pdb info into pml header
    def pml_header(self,label='Wenlin'):
        head=self.info
        head='cmd.read_pdbstr("""%s""","%s")\n' % (head.replace('\n','\\\n'),label)
        cmd="cmd.select('Chain_A','Chain A')\n"
        cmd+="cmd.hide('everything')\n"
        cmd+="cmd.show_as('cartoon','Chain_A')\n"
        cmd+="cmd.color('white','Chain_A')\n"
        cmd+="cmd.disable('Chain_A')\n"
        cmd+="cmd.center('Chain_A')\n"
        #cmd+="cmd.color('red','resi 55+44')\n"
        return head+cmd

    #pair_list is not used yet, return all_atom_dict, Calpha_dict
    def get_coordinates(self,ATOM_only=False,full_atom=False,no_water=True,Cb=False):
        """
        all_dict like dict[1]=[1.2,2.3,0.43]
        if full_dict=True, return dict[1]={'N':[1.3,2.4,2.4],atom_name:coordinate]}
        if Cb=True, return Cb distance as well. if Gly, returned Cb is None
        """
        adict={}
        bdict={}#the beta distance
        all_dict={}
        full_dict={}
        lines=self.info.split('\n')
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if ATOM_only and (not line.startswith('ATOM')):
                    continue;
                elif no_water and line[17:20].strip()=='HOH':
                    continue;
                else:
                    x=float(line[30:38])
                    y=float(line[38:46])
                    z=float(line[46:54])
                    num=int(line[22:26].strip())

                    atom_name=line[12:16].strip()
                    if full_atom:#True
                        full_dict=put_in_dict(num,full_dict,{atom_name:[x,y,z]})
                    else:
                        all_dict=put_in_dict(num,all_dict,[x,y,z])
                        if atom_name=='CA':
                            adict=put_in_dict(num,adict,[x,y,z])
                        if Cb==True:#if we need to output Cb distance
                            resi_name=line[17:20]
                            if resi_name=="GLY" and atom_name=='CA':
                                bdict=put_in_dict(num,bdict,[x,y,z])
                            elif atom_name=='CB':
                                bdict=put_in_dict(num,bdict,[x,y,z])

        if full_atom:
            return full_dict
        else:
            if Cb==False:
                return all_dict,adict
            else:#Cb==True
                return all_dict,adict,bdict


    def move(self,dist):
        "move all the atom to certain distance"
        dist=float(dist)
        new=[]
        for line in self.info.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x=float(line[30:38])+dist
                y=float(line[38:46])+dist
                z=float(line[46:54])+dist
                #print '{:>8}'.format('.3f' % x), '{:>8}'.format('.3f' % y), '{:>8}'.format('.3f' % z)
                replaced_crd='%s%s%s' % ('{:>8}'.format('%.3f' % x), '{:>8}'.format('%.3f' % y), '{:>8}'.format('%.3f' % z))
                new.append('%s%s%s' % (line[:30],replaced_crd,line[54:]))
            else:
                new.append(line)
        return '\n'.join(new)

    def cut(self,resi=None,chain='A'):
        """
        resi is the number beginning at 1 and in line[22:26]
        """
        if resi==None:
            return ''
        else:
            new=[]
            a,b=int(resi[0]),int(resi[-1])
            for line in self.info.split('\n'):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    num=int(line[22:26])
                    if num>=a and num<=b:
                        new.append(line)
            return '\n'.join(new)

    def calculate_distance(self,scale=9998,spacing=5,sort=True,as_constraint=False,ATOM_only=True,Cb=True):
        """
        if distance < scale A, take into consideration
        sort indicate sort it by  distance
        spacing==residues between two interacting residues
        scale here is now being changed to Cb (originally Ca) distance
        as_constraint will make indexs starting with 0, thus
        the output format can be used as DCA constraint format
        If Cb=True, calculate Cb distance. if False, calculate Ca
        """
        #distance cutoff in calculating closest distance
        alpha_dist_cutoff=float(scale)#A
        spacing=int(spacing)

        #compute dict
        #CA_dist righ now is not refer to CA_dist, it refered to Cb distance if Cb==True
        if Cb==False:
            atom_dict,CA_dict=self.get_coordinates(ATOM_only=ATOM_only,Cb=Cb)
        else:
            atom_dict,sth,CA_dict=self.get_coordinates(ATOM_only=ATOM_only,Cb=Cb)

        info_list=[]
        #print atom_dict['61'],atom_dict['61']
        alist=list(CA_dict.keys())
        alist.sort()
        #print alist
        blist=[]
        for i,A in enumerate(alist):
            for B in alist[i+1:]:
                #A=str(A)
                #B=str(B)
                CA_A=CA_dict[A]
                CA_B=CA_dict[B]
                CA_dist=shortest_distance(CA_A,CA_B)
                if CA_dist<alpha_dist_cutoff and abs(int(A)-int(B))>spacing:
                    blist+=[CA_dist]
                    co_A=atom_dict[A]
                    co_B=atom_dict[B]
                    co_dist=shortest_distance(co_A,co_B)
                    if as_constraint==True:
                        AA=int(A)-1
                        BB=int(B)-1
                    else:
                        AA=A
                        BB=B
                    info_list+=['%s\t%s\t%s\t%s' % (CA_dist,co_dist,AA,BB)]
        if sort:
            info_list=cmn.sort_key(blist,info_list)
        return info_list

    def interacting_pairs(self,ATOM_only=True,spacing=5,rule='CB'):
        """
        return the interacting pairs. default rule(CB) is 8A between beta residues.
        """
        pairs=[]
        if rule=='CB':
            atom_dict,sth,CB_dict=self.get_coordinates(ATOM_only=ATOM_only,Cb=True)
            alist=list(CB_dict.keys())
            for i,A in enumerate(alist):
                for B in alist[i+1:]:
                    #A=str(A)
                    #B=str(B)
                    CB_A=CB_dict[A]
                    CB_B=CB_dict[B]
                    CB_dist=shortest_distance(CB_A,CB_B)
                    if CB_dist<8 and abs(int(A)-int(B))>spacing:
                        pairs.append(tuple([int(A),int(B)]))
        return pairs

    def check_out(self):
        return self.info

    def interact_map(self):
        #2D interaction map
        alist=self.calculate_distance(spacing=0,sort=False)
        itr_ca_dict={}
        itr_closest_dict={}
        length=0
        for line in alist:
            items=line.split()
            i=int(items[-2])-1
            j=int(items[-1])-1
            length=max(i,j,length)
            itr_ca_dict['%s %s' % (i,j)]=float(items[0])
            itr_closest_dict['%s %s' % (i,j)]=float(items[1])
            #to make it symmetic
            itr_ca_dict['%s %s' % (j,i)]=float(items[0])
            itr_closest_dict['%s %s' % (j,i)]=float(items[1])
        CA_map=itmap(title='C alpha distance',length=length,interaction=itr_ca_dict,threshold=50.7758521)
        closest_map=itmap(title='C alpha distance',length=length,interaction=itr_closest_dict,threshold=50.7758521)
        return [CA_map.checkout(),closest_map.checkout()]

    def pdb2fasta(self,chain='A',ATOM_only=True):
        """
        transfer sequence in chain A to fasta format
        only consider ATOM tag
        """
        no_water=True
        resi_dict=residue_three2one_dict()
        seq=[]
        resi_numb=-999
        for line in self.info.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if ATOM_only and (not line.startswith('ATOM')):
                    continue;
                elif no_water and line[17:20].strip()=='HOH':
                    continue;
                elif line[21]==chain:
                    now_numb=int(line[22:26])
                    if now_numb!=resi_numb:
                        seq.append(resi_dict[line[17:20].strip()])
                        resi_numb=now_numb
        return ''.join(seq)

    def pdb2constraints(self,chain='A',local_range=5,move=0,ATOM_only=True,Ca_dist=None,f_horiz=None):
        """
        Calculate the distance and angle constraints between local_range
        If Ca_dist!=None, also include Ca distance constraint within Ca_dist in dist_dict
        If dssp file indicated, then only take local constraints within one element
        indexing starts with 1 (not 0)
        distant constraints between O-O,N-N,CA-CA,CB-CB,O-N,CA-O
        dist_dict like dict[(1,3)]={'O-O':float,'O-N':float}
        angle constraints for phi, psi angle
        angle is 'c(-1)-n-ca-c' and 'n(-1)-ca(-1)-c(-1)-n'
        angle_dict like dict['2phi']=float(angle)
        #move is used to shift output indexing
        """
        def parse_dist(ii,jj,move,dist_dict):
            if True:
                if True:
                    key=tuple([ii+move,jj+move])
                    if key not in dist_dict:
                        dist_dict[key]={}
                    try:
                        dist_dict[key]=put_in_dict('O-O',dist_dict[key],distance(full_dict[ii]['O'],full_dict[jj]['O']))
                    except:
                        print('weird atom in pair O%s -- O%s' % (ii,jj), file=sys.stderr)

                    try:
                        dist_dict[key]=put_in_dict('N-N',dist_dict[key],distance(full_dict[ii]['N'],full_dict[jj]['N']))
                    except:
                        print('weird atom in pair N%s -- N%s' % (ii,jj), file=sys.stderr)

                    try:
                        dist_dict[key]=put_in_dict('CA-CA',dist_dict[key],distance(full_dict[ii]['CA'],full_dict[jj]['CA']))
                    except:
                        print('weird atom in pair CA%s -- CA%s' % (ii,jj), file=sys.stderr)

                    try:#Gly
                        dist_dict[key]=put_in_dict('CB-CB',dist_dict[key],distance(full_dict[ii]['CB'],full_dict[jj]['CB']))
                    except:
                        pass

                    try:
                        dist_dict[key]=put_in_dict('O-N',dist_dict[key],distance(full_dict[ii]['O'],full_dict[jj]['N']))
                    except:
                        print('weird atom in pair O%s -- N%s' % (ii,jj), file=sys.stderr)

                    try:
                        dist_dict[key]=put_in_dict('CA-O',dist_dict[key],distance(full_dict[ii]['CA'],full_dict[jj]['O']))
                    except:
                        print('weird atom in pair CA%s -- O%s' % (ii,jj), file=sys.stderr)
            return dist_dict
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def parse_angle(ii,move,angle_dict):
            if True:
                try:
                    angle_dict=put_in_dict('%sphi' % (ii+move), angle_dict, dan(np.array(full_dict[ii-1]['C']),np.array(full_dict[ii]['N']),np.array(full_dict[ii]['CA']),np.array(full_dict[ii]['C'])))
                except:
                    print('weird atom for phi angle %s' % ii, file=sys.stderr)

                try:
                    angle_dict=put_in_dict('%spsi' % (ii+move), angle_dict, dan(np.array(full_dict[ii-1]['N']),np.array(full_dict[ii-1]['CA']),np.array(full_dict[ii-1]['C']),np.array(full_dict[ii]['N'])))
                except:
                    print('weird atom for psi angle %s' % ii, file=sys.stderr)
            return angle_dict
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        angle_dict,dist_dict={},{}
        full_dict=self.get_coordinates(ATOM_only=ATOM_only,full_atom=True)
        k=max(full_dict.keys())
        left=min(full_dict.keys())#usually, left is 1
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if f_horiz!=None:#if horiz file indicated, only take the range of s. element
            regions=hz(f_horiz).ss_region(combine=True)
            for region in regions:
                a,b=region
                for i in range(a,b):
                    ii=i+left
                    if ii>1:
                        angle_dict=parse_angle(ii,move,angle_dict)
                    for j in range(i+1,min((i+1+local_range),b)):
                        jj=j+left
                        if jj<=k:
                            dist_dict=parse_dist(ii,jj,move,dist_dict)
        else:#if no s.s. indicated
            for i in range(k):
                ii=i+left
                if ii>1:
                    angle_dict=parse_angle(ii,move,angle_dict)
                for j in range(i+1,i+1+local_range):
                    jj=j+left
                    if jj<=k:
                        dist_dict=parse_dist(ii,jj,move,dist_dict)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if Ca_dist!=None:#if remote distance constraint needed, add here
            for i in range(k):
                ii=i+left
                for j in range(i+1+local_range,k):#jump 1 residue
                    jj=j+left
                    raw_dist=distance(full_dict[ii]['CA'],full_dict[jj]['CA'])
                    if raw_dist<=Ca_dist:
                        key=tuple([ii+move,jj+move])
                        if key not in dist_dict:
                            dist_dict[key]={}
                        dist_dict[key]=put_in_dict('CA-CA',dist_dict[key],raw_dist)
        return dist_dict,angle_dict


    def pdb2noe(self,dist_weight=3,chain='A',local_range=5,move=0,ATOM_only=True,Ca_dist=None,f_horiz=None):
        """
        take all the constraints in pdb file and make them as noe info
        if Ca_dist indicated, take remote constraints between Ca_dist
        right now, only designed to take real constraint from structure
        return lines (dist,angle)
        """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #used in Rosetta averaging
        #def extracting(alist,digit=2):
        #    b=np.array(alist)
        #    return round(b.mean(),digit),round(b.min(),digit),round(b.std(),digit)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        resi_dict=residue_one2three_dict()
        seq=self.pdb2fasta(chain,ATOM_only)
        dist,angle=[],[]
        dist_dict,angle_dict=self.pdb2constraints(chain,local_range,move,ATOM_only=ATOM_only,Ca_dist=Ca_dist,f_horiz=f_horiz)
        for resi_pair in dist_dict:
            #6 keys
            i,j=resi_pair
            for atom_pair in dist_dict[resi_pair]:
                data=dist_dict[resi_pair][atom_pair]
                #Mean,Min,SD=extracting(data)
                Mean=round(data[0],1)
                Min=Mean-1
                SD='1.0'
                atoms=atom_pair.split('-')
                dist.append('assign (resid %s and name %s) (resid %s and name %s)  %s %s %s weight %s ! %s %s' %
                            (i,atoms[0],j,atoms[-1],Mean,Min,SD,dist_weight,resi_dict[seq[i-1]],resi_dict[seq[j-1]]))

        for key in angle_dict:
            i=int(key.strip('phsi'))
            data=angle_dict[key]
            #Mean,Min,SD=extracting(data)
            Mean=round(data[0],1)
            SD='20.0'
            #angle is phi'c(-1)-n-ca-c' and psi'n(-1)-ca(-1)-c(-1)-n'
            if 'phi' in key:
                angle.append('assign (resid %s and name c) (resid %s and name n) (resid %s and name ca) (resid %s and name c)  1.0 %s %s 2' %
                    (i-1,i,i,i,Mean,SD))
            elif 'psi' in key:
                angle.append('assign (resid %s and name n) (resid %s and name ca) (resid %s and name c) (resid %s and name n)  1.0 %s %s 2' %
                    (i-1,i-1,i-1,i,Mean,SD))
        return dist,angle

    def pair2pml(self,pairs,colored_pairs=20, selected_pairs=50, label='STH'):
        """
        map interaction pairs to structure, generate pml file.
        pairs should be like [[1,2],[3,4]]
        """
        info=self.pml_header(label=label)
        #use the old function above
        info+=color_top(pairs,colored_pairs)
        info+=select_top(pairs,selected_pairs)
        return info

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#x,y is a list of three element
def distance(x,y):
    sq=[(x[i]-y[i])*(x[i]-y[i]) for i in range(3)]
    return math.sqrt(sum(sq))

