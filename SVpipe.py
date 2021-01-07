#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import SVread as read
import SVsee as see
import argparse
import time
import sys

def _regfilt(start,stop,reg):
    if reg[0]<=start<=reg[1] and reg[0]<=stop<=reg[1] or start<=reg[0]<=stop and start<=reg[1]<=stop:
       return True

def _limfilt(rvars,**kwargs):
    #removes variants that lie outside of bounds from rvars
    reg=kwargs.get('reg',None)
    llim=kwargs.get('llim',None)
    ulim=kwargs.get('ulim',None)
    newvars=[]
    if reg:
       for x in rvars:
           if read._lencheck(x.length,llim=llim,ulim=ulim) and _regfilt(x.start,x.stop,reg):
              newvars.append(x)
    else:
       for x in rvars:
           if read._lencheck(x.length,llim=llim,ulim=ulim):
              newvars.append(x)
    return newvars

def _vdraw2(vlist,Plot):
    #modified version of _vdraw from see, takes variant list immediately
    vdict=see._stack(vlist,0.015)
    max_h=max(vdict.keys())+0.015
    cdict=see._cdict(vlist,intype='vcf')
    #max_h gives lower y-boundary of highest row!
    Plot.ax1.set_ylim(0,max_h)
    see._fbox(vdict,Plot.ax1,0.0135,cdict)

def _annodraw2(annolist,format,Plot):
    #modified version of see's _annodraw, takes feature list immediately
    # -> no tupperware, relevant for GUI
    if annolist:
       annodict=see._stack(annolist,0.015)
       max_h=max(annodict.keys())+0.015
       Plot.ax2.set_ylim(0,max_h)
       if format=='gff':
          cdict=see._cdict(annolist,intype=format)
          see._fbox(annodict,Plot.ax2,0.0135,cdict)
       else:
          for key,value in annodict.items():
              for subval in value:
                  p=patches.Rectangle((subval.start, key), subval.length, 0.0135, color='b')
                  Plot.ax2.add_patch(p)

if __name__=='__main__':
   start_time = time.time()

   parser=argparse.ArgumentParser(description='Full structural variant visualization from raw VCF file to figure')
   parser.add_argument('--vcf',help='Provide VCF file')
   parser.add_argument('-C','--contigname',type=str,nargs=1,help='Provide contig/chromosome ID by which to filter SV records; multiple IDs given space-separated; use only one ID for visualization')
   parser.add_argument('-T','--type',type=str,nargs='+',help="Provide type(s) by which to filter SV records i.e. 'DEL'; multiple types given space-separated; 'BND' for complex rearrangments with breakends")
   parser.add_argument('-S','--conttype',type=str,nargs='+',help="Verbose type filter: types entered will be excluded. Same usage as -T")
   parser.add_argument('--ulim',type=int,help="Provide upper length limit for SVs [bp] (inclusive)")
   parser.add_argument('--llim',type=int,help="Provide lower length limit for SVs [bp] (inclusive)")
   parser.add_argument('--gff',dest='mod',help='Provide contig annotation as GFF file; must have .gff or .GFF extension')
   parser.add_argument('--bed',dest='mod',help='Provide contig annotation as BED file; must have .bed or .BED extension')
   parser.add_argument('--blast',dest='mod',help='Provide contig annotation as BLAST output file (outfmt 6)')
   parser.add_argument('-R','--anno_type',type=str,nargs='+',help='Provide type(s) of annotations to be included (only applicable on GFF annotation files)')
   parser.add_argument('-RC','--anno_cont',type=str,nargs='+',help='Complementary type filter...analogous to -S (only applicable on GFF annotation files)')
   parser.add_argument('-RL','--anno_lower',type=int,help='Provide lower length limit for annotation features (in bp)')
   parser.add_argument('-RU','--anno_upper',type=int,help='Provide upper length limit for annotation features (in bp)')    
   parser.add_argument('--heat',type=int,help='Generate heatmap with given window size (in bp)')
   parser.add_argument('--reg',nargs=2,type=int,help='Specify region bounds for more zoomed-in look; bounds given whitespace-separated i.e. 1000 1500')
   parser.add_argument('--out',type=str,help='Provide desired name of graphic output')
   parser.add_argument('--outfmt',choices=['pdf','svg','png'],default='png',help='Select output format, default=png')
   parser.add_argument('--ov', action='store_true', help='Choose to display an overview of SVs in the provided VCF file; --ov only will print assembly-level information; using --C option will print single-contig-level information')
   args=parser.parse_args()

   #read vcf file into object list, filter by various filtering options
   variants=read._readfile(args.vcf, cfil=args.contigname, tfil=args.type, vfil=args.conttype)
   variants.sort(key=lambda x: x.start)
   nvars=0
   if args.llim or args.ulim:
      nvars=_limfilt(variants, llim=args.llim, ulim=args.ulim)
   else:
      nvars=variants

   if args.ov:
      if args.contigname and nvars:
         read._getstats(nvars)
      elif args.contigname:
         read._getstats(variants)
      else:
         read._contig_overview(args.vcf)
      sys.exit()

   #plotting
   cname=variants[0].cont
   clen=see._chrxlen(args.vcf,cname)
   plt.autoscale(False)
   if args.mod and args.heat:
      plot = see.Plot(cname,clen,reg=args.reg, anno=args.mod, heat=args.heat)
      _vdraw2(nvars,plot)
      see._annodraw(args.mod,cname,plot,reg=args.reg,a_tfil=args.anno_type,a_lo=args.anno_lower,a_up=args.anno_upper,a_co=args.anno_cont)
      see._heatmap(nvars,clen,args.heat,plot,reg=args.reg)
   elif args.mod:
      plot = see.Plot(cname,clen,reg=args.reg, anno=args.mod, heat=args.heat)
      _vdraw2(nvars,plot)
      see._annodraw(args.mod,cname,plot,reg=args.reg,a_tfil=args.anno_type,a_lo=args.anno_lower,a_up=args.anno_upper,a_co=args.anno_cont)
   elif args.heat:
      plot = see.Plot(cname,clen,reg=args.reg, anno=args.mod, heat=args.heat)
      _vdraw2(nvars,plot)
      see._heatmap(nvars,clen,args.heat,plot,reg=args.reg)
   else:
      plot = see.Plot(cname,clen,reg=args.reg, anno=args.mod, heat=args.heat)
      _vdraw2(nvars,plot)

   #plt.tight_layout(pad=1,h_pad=0.35)

   plt.savefig(args.out,format=args.outfmt)
   print('-- Run time: {:.2f} --'.format(time.time()-start_time))
