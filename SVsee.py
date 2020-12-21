#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib._color_data as mcd
from matplotlib import gridspec
from matplotlib import ticker
from matplotlib import style

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from collections import OrderedDict
import argparse
import time
import re
import random

class Variant:
#bundles all relevant information per vcf file line
    def __init__(self,cont,start,length,type):
        self.cont=cont
        self.start=int(start)
        self.length=int(length)
        self.end=self.start+self.length
        self.type=type

    def _add2list(self,target_list):
        target_list.append(self)

class Feature:
#bundles information from gff file line
    def __init__(self,cont,type,start,stop,info):
        self.cont=cont
        self.type=type
        self.start=int(start)
        self.end=int(stop)
        self.length=self.end-self.start
        self.info=info

    def _add2list(self,target_list):
        target_list.append(self)

##implement more input file formats##
class Bed:
    def __init__(self,chrom,start,end):
        self.cont=chrom
        self.start=int(start)
        self.end=int(end)
        self.length=self.end-self.start

    def _add2list(self,targetlist):
        targetlist.append(self)

##blast outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore##
class Blasthit:
    def __init__(self,qseqid,sseqid,qstart,qend,sstart,send):
        self.qseqid=qseqid
        self.qstart=int(qstart)
        self.qend=int(qend)
        self.cont=sseqid
        self.start=int(sstart)
        self.end=int(send)
        self.length=abs(self.end-self.start)

    def _add2list(self,targetlist):
        targetlist.append(self)


##########plot-setup########
def _chrx(infile):
    with open(infile,'r') as inhandle:
         line=inhandle.readline()
    return line.split('\t')[0]
    inhandle.close()

def _chrxlen(vcf,cname):
    for line in open(vcf,'r'):
        if line[0]=='#' and cname in line:
           lline=line
           break
    clen=re.findall(',length=[0-9]+',lline)
    clen=int(clen[0].split('=')[1])
    return clen

class Plot:
    def __init__(self,cname,clen,**kwargs):
        self.reg=kwargs.get('reg',None)
        self.anno=kwargs.get('anno',None)
        self.heat=kwargs.get('heat',None)
        self.cname=cname
        self.clen=clen
        self._setupplot()

    def _setupplot(self):
        plt.rcParams["figure.figsize"] = 15,7
        style.use('ggplot')
        self.fig = plt.figure()
        if self.anno and self.heat:
           spec = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=[7,7,1,1])
           self.ax1 = self.fig.add_subplot(spec[0])
           self.ax2 = self.fig.add_subplot(spec[1],sharex=self.ax1)
           self.ax3 = self.fig.add_subplot(spec[2],sharex=self.ax1)
           self.ax4 = self.fig.add_subplot(spec[3],sharex=self.ax1)
           self.ax3.grid(False)
           self.ax4.set_xlabel('Nucleotide position [b]')
           self.num_ax=4
        elif self.anno:
           spec = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[8,8])
           self.ax1 = self.fig.add_subplot(spec[0])
           self.ax2 = self.fig.add_subplot(spec[1],sharex=self.ax1)
           self.ax2.set_xlabel('Nucleotide position [b]')
           self.num_ax=2
        elif self.heat:
           spec = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[14,1,1],)
                                    #left=0.1,
                                    #right=0.9)
           self.ax1 = self.fig.add_subplot(spec[0])
           self.ax2 = self.fig.add_subplot(spec[1],sharex=self.ax1)
           self.ax2.grid(False)
           self.ax3 = self.fig.add_subplot(spec[2],sharex=self.ax1)
           self.ax3.set_xlabel('Nucleotide position [b]')
           self.num_ax=3
        else:
           self.ax1 = self.fig.add_subplot(1,1,1)
           self.ax1.set_xlabel('Nucleotide position [b]')
           self.ax2=None
           self.num_ax=1
        self._retxlims()
        self._setspineslabels()
        self.ax1.set_title(f'SV Distribution on {self.cname}',fontsize=17)
        
        plt.tight_layout(pad=3.0,h_pad=0.35)

    def _retxlims(self):
        if self.reg:
           lo,hi=self.reg[0],self.reg[1]
        else:
           lo,hi=0,self.clen

        if self.ax2:
           self.ax1.set_xlim(lo,hi)
           self.ax2.set_xlim(lo,hi)
        else:
           self.ax1.set_xlim(lo,hi)

    def _setspineslabels(self):
        self.ax1.spines['top'].set_visible(False)  #equivalent to ax1.set_frame_on(False)
        self.ax1.spines['right'].set_visible(False)
        self.ax1.spines['left'].set_visible(False)

        self.ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        self.ax1.axes.get_yaxis().set_visible(False)

        if self.ax2:
           self.ax2.spines['top'].set_visible(False)  #equivalent to ax1.set_frame_on(False)
           self.ax2.spines['right'].set_visible(False)
           self.ax2.spines['left'].set_visible(False)
           self.ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
           self.ax2.axes.get_yaxis().set_visible(False)

#########data-wrangling############
def _cfit(line,cname):
    if cname==None:
       return True
    else:
       return(True if line[0:len(cname)]==cname else False)

def _tupperware(input,cname,**kwargs):
    type=kwargs.get('type',None)
    reg=kwargs.get('reg',None)
    if type=='gff' or type=='GFF':
       indices=[0,2,3,4,8]
       getclass=Feature
    elif type=='bed' or type=='BED':
       indices=[0,1,2]
       getclass=Bed
    elif type=='blast':
       indices=[0,1,6,7,8,9]
       getclass=Blasthit
    elif type=='variants':
       indices=[0,1,2,3]
       getclass=Variant
    else:
       indices=[0,1,6,7,8,9]
       getclass=Blasthit
    targetlist=[]
    for line in open(input,'r'):
        if line[0]!='#':
           a=getclass(*[line.split('\t')[x] for x in indices])
           if a.cont==cname:
              if reg:
                 if reg[0]<=a.start<=reg[1]:
                    a._add2list(targetlist)
              else:
                 a._add2list(targetlist)
    return targetlist

####change cdict to avoid random color assignment#########################
def _cdict(features,**kwargs):
    #takes the primary object list of variants or features!
    in_type=kwargs.get('intype',None)
    cdict=OrderedDict()
    cols = list(mcd.XKCD_COLORS.values())
    if in_type=='vcf':
       full_dict={'DEL':cols[15],'INS':cols[6],'DUP':'black',
              'INV':cols[18],'BND':'red','CNV':cols[33],
              'DUP:TANDEM':cols[83],'DEL:ME':cols[93],'INS:ME':cols[3]}
       types=set([x.type for x in features])
    elif in_type=='gff':
        full_dict={'region':cols[15],'gene':cols[6],'exon':'black',
              'CDS':cols[18],'snRNA':'red','mRNA':cols[33],
              'transcript':cols[83],'cDNA_match':cols[93],'snoRNA':cols[3],
              'lnc_RNA':cols[222],'pseudogene':cols[211],
              'primary_transcript':cols[233],'tRNA':cols[196],'miRNA':cols[133]}
        types=set([x.type for x in features])
    for type in types:
        if type in full_dict.keys():
           cdict[type]=full_dict[type]
        elif type not in full_dict.keys():
             cdict[type]=random.choice(cols)
    return cdict

def _legend(ax,xdict):
    patchList = []
    for key in xdict:
        data_key = patches.Patch(color=xdict[key], label=key)
        patchList.append(data_key)
    ax.legend(handles=patchList, loc='upper left', bbox_to_anchor=(0,1.0), ncol=5,fancybox=True)

def _fbox(hdict,ax,height,cdict):
    for key,value in hdict.items():
        for subval in value:
            if 'info' in dir(subval):
               p=patches.Rectangle((subval.start, key),
                                subval.length, height,
                                color=cdict[subval.type],
                                fill=True, picker=5, label=subval.info)
            else:
               p=patches.Rectangle((subval.start, key),
                                subval.length, height,
                                color=cdict[subval.type],
                                fill=True, picker=5)
            ax.add_patch(p)
    _legend(ax,cdict)

def _laststop(hdict):
    return sorted([(x,y[-1].end) for x,y in hdict.items()],key=lambda x: x[1])

def _stack(vlist,step):
#works with both variant and feature-type objects!
    hdict={}
    ini_height=0
    hdict[ini_height]=[vlist[0]]
    for i in range(1,len(vlist)):
        for x in _laststop(hdict):
            if vlist[i].start>x[1]:
               hdict[x[0]]+=[vlist[i]]
            else:
               hdict[ini_height+step]=[vlist[i]]
               ini_height+=step
               break
            break
    return hdict

def _vdraw(infile,cname,Plot,**kwargs):
    reg=kwargs.get('reg',None)
    ret=kwargs.get('ret',None)
    vlist=_tupperware(infile,cname,type='variants',reg=reg)
    cdict=_cdict(vlist,intype='vcf')
    vdict=_stack(vlist,0.015)
    max_h=max(vdict.keys())+0.015
    Plot.ax1.set_ylim(0,max_h)
    _fbox(vdict,Plot.ax1,0.0135,cdict)
    if ret:
       return vlist

def _annodraw(input,cname,Plot,**kwargs):
    reg=kwargs.get('reg',None)
    tfil=kwargs.get('a_tfil',None)
    a_lo=kwargs.get('a_lo',None)
    a_up=kwargs.get('a_up',None)
    a_co=kwargs.get('a_co',None)
    def _filterbylength(input,lower_limit,upper_limit):
        output=[]
        if lower_limit and upper_limit:
           for x in input:
               if lower_limit<=x.length<=upper_limit:
                  output.append(x)
        elif lower_limit:
           for x in input:
               if lower_limit<=x.length:
                  output.append(x)
        elif upper_limit:
           for x in input:
               if upper_limit>=x.length:
                  output.append(x)
        else:
           output=input[:]
        return output   
           
    ext=re.split('\.',input)[-1]
    anno_list=_tupperware(input,cname,type=ext,reg=reg)
    anno_list=_filterbylength(anno_list, a_lo, a_up)
    if ext=='gff' and tfil or ext=='GFF' and tfil:
       new_anno=[x for x in anno_list if x.type in tfil]
       anno_list=sorted(new_anno,key=lambda x: x.start)
    elif ext=='gff' and a_co or ext=='GFF' and a_co:
       new_anno=[x for x in anno_list if x.type not in a_co]
       anno_list=sorted(new_anno,key=lambda x: x.start)
    
    anno_dict=_stack(anno_list,0.015)
    max_h=max(anno_dict.keys())+0.015
    Plot.ax2.set_ylim(0,max_h)
    if ext=='gff' or ext=='GFF':
       cdict=_cdict(anno_list,intype=ext)
       _fbox(anno_dict,Plot.ax2,0.0135,cdict)
    else:
       for key,value in anno_dict.items():
           for subval in value:
               p=patches.Rectangle((subval.start, key),subval.length, 0.0135)
               Plot.ax2.add_patch(p)

#####heatmap###############
def _ranges(edict,lo,hi,step):
    for i in range(lo+step,hi+step,step):
        edict[i]=0

def __roundup(x,step):
    return x if x % step == 0 else x + step - x % step

def _altrangehits(elem,ranges,lo,hi,step):
    start=__roundup(elem.start,step)
    stop=__roundup(elem.end,step)
    for i in range(start,stop+step,step):
        try:
           ranges[i]+=1
        except KeyError:
           pass

def _arrpoints(lo,hi,step):
    lo=int(lo+(step/2))
    hi=int(hi+(step/2))
    return [x for x in range(lo,hi,step)]

def _heat(features,lo,hi,step):
    ranges=OrderedDict([])
    _ranges(ranges,lo,hi,step)
    [_altrangehits(x,ranges,lo,hi,step) for x in features]
    return ranges

def _heatmap(features,clen,step,Plot,**kwargs):
    reg=kwargs.get('reg',None)
    if reg:
       lo=reg[0]
       hi=reg[1]
    else:
       lo=0
       hi=clen

    ranges=_heat(features,lo,hi,step)

    import numpy as np
    x=np.array(_arrpoints(lo,hi,step))
    y=np.array([value for value in ranges.values()])
    maxy=max([value for value in ranges.values()])
    medy=__roundup((maxy/2),1)
    if Plot.num_ax==4:
       ax_heat=Plot.ax3
       ax_curve=Plot.ax4
    elif Plot.num_ax==3:
       ax_heat=Plot.ax2
       ax_curve=Plot.ax3

    extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
    ax_heat.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    _heatmap=ax_heat.imshow(y[np.newaxis,:], cmap="plasma", aspect="auto", extent=extent)
    ax_heat.set_yticks([])
    ax_heat.set_xlim(extent[0], extent[1])
 
    axins = inset_axes(ax_heat,
                           width="1%",  # width = 5% of parent_bbox width
                           height="100%",  # height : 50%
                           loc='center left',
                           bbox_to_anchor=(-0.011, 0, 1, 1),
                           bbox_transform=ax_heat.transAxes,
                           borderpad=0,
                           )
    
    cbar = ax_heat.figure.colorbar(_heatmap,
                                   cax=axins,
                                   orientation='vertical',
                                   fraction=.25,
                                   ticks=[0,maxy])
    cbar.ax.yaxis.set_ticks_position('left') 
    cbar.ax.yaxis.set_tick_params(labelsize=9)
    
    ax_curve.plot(x,y)
    ax_curve.set_xlim(extent[0], extent[1])
    ax_curve.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax_curve.set_yticks([0,medy,maxy])
    ax_curve.yaxis.set_tick_params(labelsize=9)
    ax_curve.set_ylabel(f'Variant count per\n{step} bp',
                        fontsize='small',
                        labelpad=-0.4,
                        va='bottom')


if __name__ == "__main__":
   start_time = time.time()

   parser=argparse.ArgumentParser(description='Single-chromsome visualization of SV distribution based on SVread.py output')
   parser.add_argument('--IN',help='Provide input file (SVread.py output file)')
   parser.add_argument('--vcf',help='Provide VCF file corresponding to input file')
   parser.add_argument('--gff',dest='mod',help='Provide contig annotation as gff file; must have .gff or .GFF extension')
   parser.add_argument('--bed',dest='mod',help='Provide contig annotation as bed file; must have .bed or .BED extension')
   parser.add_argument('--blast',dest='mod',help='Provide contig annotation as blast output file (outfmt 6)')
   parser.add_argument('-R','--anno_type',type=str,nargs='+',help='Provide type(s) of annotations to be included (only applicable on GFF annotation files)')
   parser.add_argument('-RC','--anno_cont',type=str,nargs='+',help='Complementary type filter, selected type(s) will be excluded, multiple arguments given space-separated')
   parser.add_argument('-RL','--anno_lower',type=int,help='Provide lower length limit for annotation features (in bp)')
   parser.add_argument('-RU','--anno_upper',type=int,help='Provide upper length limit for annotation features (in bp)')
   parser.add_argument('--heat',type=int,help='Generate heatmap with given window size (in bp)')
   parser.add_argument('--reg',nargs=2,type=int,help='Specify region bounds (in bp) for more zoomed-in look; bounds given whitespace-separated i.e. 1000 1500')
   parser.add_argument('--out',type=str,help='Provide desired name of output file')
   parser.add_argument('--outfmt',choices=['pdf','svg','png'],default='png',help='Select output format, default=png')
   args=parser.parse_args()

   cname=_chrx(args.IN)
   clen=_chrxlen(args.vcf,cname)

   if args.mod and args.heat:
      plot = Plot(cname,clen,reg=args.reg, anno=args.mod, heat=args.heat)
      variants=_vdraw(args.IN,cname,plot,reg=args.reg,ret=True)
      _annodraw(args.mod,cname,plot,reg=args.reg,a_tfil=args.anno_type,a_lo=args.anno_lower,a_up=args.anno_upper,a_co=args.anno_cont)
      _heatmap(variants,clen,args.heat,plot,reg=args.reg)
   elif args.mod:
      plot = Plot(cname,clen,reg=args.reg, anno=args.mod)
      _vdraw(args.IN,cname,plot,reg=args.reg)
      _annodraw(args.mod,cname,plot,reg=args.reg,a_tfil=args.anno_type,a_lo=args.anno_lower,a_up=args.anno_upper,a_co=args.anno_cont)
   elif args.heat:
      plot = Plot(cname,clen,reg=args.reg, heat=args.heat)
      variants=_vdraw(args.IN,cname,plot,reg=args.reg,ret=True)
      _heatmap(variants,clen,args.heat,plot,reg=args.reg)
   else:
      plot = Plot(cname,clen,reg=args.reg)
      _vdraw(args.IN,cname,plot,reg=args.reg)

   plt.savefig(args.out,format=args.outfmt)
   print('-- Run time: {:.2f} --'.format(time.time()-start_time))
