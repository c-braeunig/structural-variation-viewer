#!/usr/bin/env python3

import time
import sys
import re
import argparse
import statistics


class Variant:
    def __init__(self,cont,start,length,type,seq):
        self.type=type
        self.cont=cont.rstrip().lstrip()
        self.start=int(start)
        self.length=int(length) if length>0 else 1
        self.end=self.start+self.length
        self.seq=seq

    def add2list(self,targetlist):
        targetlist.append(self)

class Link(Variant):
    def __init__(self,Variant,target,tstart):
        self.cont=Variant.cont
        self.start=Variant.start
        self.length=Variant.length
        self.end=Variant.start+Variant.length
        self.type=Variant.type
        self.seq=Variant.seq
        self.target=target
        self.tstart=int(tstart)
        self.tend=self.tstart+Variant.length

    def add2list(self,targetlist):
        targetlist.append(self)

class Contig:
    def __init__(self,name,length):
        self.name=name
        self.label=name[:5]
        self.start=0
        self.length=length

    def add2list(self,targetlist):
        targetlist.append(self)

#############generate variant object#########
def _regex():
    global end_rex
    global size_rex
    global seq_rex
    global type_rex
    global indices
    end_rex=re.compile(';END=[0-9]+')
    size_rex=re.compile(';SVLEN=-?[0-9]+')
    seq_rex=re.compile(';SEQ=[agctAGCT]+')
    type_rex=re.compile(';SVTYPE=<?[A-Z]+>?')
    indices=[0,1,3,4,7]

def _sepcolumns(line):
    indices=[0,1,7,4]
    _split=line.split
    columns=[_split()[x] for x in indices]
    return columns

def _end(columns, rex):
    try:
       end=rex.findall(columns[4])
       end=int(end[0].split('=')[1])
    except IndexError:
       end='.'
    return end

def _size(columns,size_re,seq_re):
    try:
       size=size_re.findall(columns[4])
       size=abs(int(size[0].split('=')[1]))
       if size==999999999:
          raise IndexError
    except IndexError:
       end=_end(columns,end_rex)
       size=end-int(columns[1])
       if size==0:
          size=len(_seq(columns,seq_re))
    return size

def _seq(columns,seq_re):
    try:
       seq=seq_re.findall(columns[4])
       seq=seq[0].split('=')[1]
       return seq
    except IndexError:
       if '<' not in columns[4]:
          return(columns[2] if len(columns[2])>len(columns[3]) else columns[3])

def _type(columns,type_re):
    try:
       type=type_re.findall(columns[4])
       type=type[0].split('=')[1]
    except IndexError:
       type='.'
    return type

###########accumulate variants into a list###########
def _apconds(element,targetlist,**kwargs):
    #new version: first generate variant object per line, then see if excluded or not
    tfil=kwargs.get('tfil',None)
    cfil=kwargs.get('cfil',None)
    vfil=kwargs.get('vfil',None)
    if element[0]!='#':
       columns=[element.split()[x] for x in indices]
       a=Variant(columns[0],int(columns[1]),_size(columns,size_rex,seq_rex),_type(columns,type_rex),_seq(columns,seq_rex))
       if cfil:
          if a.cont==cfil[0]:
             if tfil:
                if a.type in tfil:
                   a.add2list(targetlist)
             elif vfil:
                if a.type not in vfil:
                   a.add2list(targetlist)
             else:
                a.add2list(targetlist)
       else:
          if tfil:
             if a.type in tfil:
                a.add2list(targetlist)
          elif vfil:
             if a.type not in vfil:
               a.add2list(targetlist)
          else:
               a.add2list(targetlist)

def _readfile(input,**kwargs):
    tfil=kwargs.get('tfil',None)
    cfil=kwargs.get('cfil',None)
    vfil=kwargs.get('vfil',None)
    variants=[]
    _regex()
    with open(input,'r') as inhandle:
         [_apconds(line,variants,cfil=cfil,tfil=tfil,vfil=vfil) for line in inhandle]
    inhandle.close()
    return variants

##########process list for flag########
def _labeldict(variants):
        vardict={}
        for elem in variants:
            if elem.type not in vardict:
               vardict[elem.type]=[elem.length]
            else:
               vardict[elem.type].append(elem.length)
        [x.sort() for x in vardict.values()]
        return vardict

def _getstats(variants):
    def _printstats(key,value):
        try:
           low=min(value)
           high=max(value)
           avg=sum(value)/len(value)
           med=statistics.median(value)
           print(f"{key}: {len(value)}  --- MIN: {low}, MAX: {high}, AVG: {avg:.1f}, MED: {med:.0f}")
        except ValueError:
           pass

    vardict=_labeldict(variants)
    [x.sort for x in vardict.values()]
    print(f'Total number of SVs: {len(variants)}')
    [_printstats(x,y) for x,y in vardict.items()]
    
def _contig_overview(vcf_file):
    contig_dict={}
    for line in open(vcf_file,'r'):
        if line[0]!='#':
           contig=line.split()[0]
           if contig not in contig_dict:
              contig_dict[contig]=1
           elif contig in contig_dict:
              contig_dict[contig]+=1
    print(f'Overview of total SVs on {len(contig_dict.keys())} contigs')
    [print(f'{contig}: {total}') for contig,total in contig_dict.items()]


############process list for list output#########
def _gettpos(v):
    return v.seq.split(':')[1].rstrip('[]N')

def _gettname(v):
    return v.seq.split(':')[0].lstrip('[]N')

def _sortedlist(list,sorter):
    if sorter=='asc_size':
       return sorted(list, key=lambda x: x.length)
    elif sorter=='des_size':
       return sorted(list, key=lambda x: x.length, reverse=True)
    elif sorter=='type':
       return sorted(list, key=lambda x: x.type)
    elif sorter=='type_asc_size':
       return sorted(sorted(list, key=lambda x: x.length), key=lambda x: x.type)
    elif sorter=='type_des_size':
       return sorted(sorted(list, key=lambda x: x.length, reverse=True), key=lambda x: x.type)
    elif sorter=='CR':
       return sorted(sorted(list, key=_gettpos), key=_gettname)
    elif not sorter:
       return sorted(sorted(list, key=lambda x: x.start), key=lambda x: x.cont)

###########process list for circos################
class Circos:
    def __init__(self,variants,vcf):
        self.variants=variants
        self.vcf=vcf
        self.links=self._links()
        self.karyotype=self._karyotype()
        self._karyo_out()
        self._links_out()

    def _links(self):
        links=[]
        for x in self.variants:
            if x.type=='BND':
               y=Link(x,_gettname(x),_gettpos(x))
               y.add2list(links)
        return links

    def _karyotype(self):
        targets=set([x.target for x in self.links])
        targetlist=[]
        lrex=re.compile(',length=[0-9]+')
        for line in open(self.vcf,'r'):
            for target in targets:
                if line[0]=='#' and target in line:
                   lline=line
                   tlen=int(lrex.findall(lline)[0].split('=')[1])
                   a=Contig(target,tlen)
                   a.add2list(targetlist)
                   break
        return targetlist

    def _karyo_out(self):
        inhandle=self.vcf.rstrip('.vcf')+'.karyo'
        with open(inhandle,'w') as karyo:
             for x in self.karyotype:
                 karyo.write(f'chr - {x.name} {x.label} {x.start} {x.length} black\n')
        karyo.close()

    def _links_out(self):
        inhandle=self.vcf.rstrip('.vcf')+'.links'
        with open(inhandle,'w') as linko:
             for x in self.links:
                 linko.write(f'{x.cont} {x.start} {x.end} {x.target} {x.tstart} {x.tend}\n')
        linko.close()

############output################################
def _lencheck(length,**kwargs):
    ulim=kwargs.get('ulim',None)
    llim=kwargs.get('llim',None)
    if llim and ulim:
       return(True if llim<=length<=ulim else False)
    elif llim:
       return(True if llim<=length else False)
    elif ulim:
       return(True if length<=ulim else False)
    else:
       return True

def _listputter(list,**kwargs):
    ulim=kwargs.get('ulim',None)
    llim=kwargs.get('llim',None)
    nos=kwargs.get('nosize',None)
    if nos:
       [print(f'{x.cont}\t{x.start}\t{x.type}\t{x.seq}') for x in list if _lencheck(x.length,llim=llim,ulim=ulim)]
    else:
       [print(f'{x.cont}\t{x.start}\t{x.length}\t{x.type}\t{x.seq}') for x in list if _lencheck(x.length,llim=llim,ulim=ulim)]

def _list2file(list,output,**kwargs):
    ulim=kwargs.get('ulim',None)
    llim=kwargs.get('llim',None)
    nos=kwargs.get('nosize',None)
    with open(output,'w') as outhandle:
         if nos:
            [outhandle.write(f'{x.cont}\t{x.start}\t{x.type}\t{x.seq}\n') for x in list if _lencheck(x.length,llim=llim,ulim=ulim)]
         else:
            [outhandle.write(f'{x.cont}\t{x.start}\t{x.length}\t{x.type}\t{x.seq}\n') for x in list if _lencheck(x.length,llim=llim,ulim=ulim)]
    outhandle.close()

def _seq2file(list,output,**kwargs):
    ulim=kwargs.get('ulim',None)
    llim=kwargs.get('llim',None)
    with open(output,'w') as outhandle:
         [outhandle.write(f'>{x.cont}:{x.start},{x.length},{x.type}\n{x.seq}\n') for x in list if '<' not in x.seq and _lencheck(x.length,llim=llim,ulim=ulim)]
    outhandle.close()

if __name__ == "__main__":
   start_time = time.time()

   parser=argparse.ArgumentParser(description='Parse and derive SV information from .vcf output of variant calling software; --ov generates a superficial overview of the listed SV types and their number of occurrences. Using --sort, --C, --T, --S, --llim and --ulim, the SVs of intereste can be isolated and printed to stdout or to a file ( --out) in a bed-like file format. For further information on the vcf format consult https://samtools.github.io/hts-specs/VCFv4.2.pdf')
   parser.add_argument('--vcf',help='Provide variant call file in VCF format')
   parser.add_argument('--sort',choices=['asc_size','des_size','type','type_asc_size','type_des_size','CR'],help='Provide sorting criterium to organize SVs; SVs sorted by position on chromosome/contig by default; no sorting recommended for downstream SVsee use; Use "CR" when filtering by "BND" to sort target locations')
   parser.add_argument('--ov',action='store_true',help='Choose to display an overview of SVs in the provided vcf file; --ov only will print assembly-level information; using --C option will print single-contig-level information')
   parser.add_argument('-C','--contigname',type=str,nargs=1,help='Provide contig/chromosome ID by which to filter SV records; multiple IDs given space-separated; use only one for visualization')
   parser.add_argument('-T','--type',type=str,nargs='+',help="Provide type(s) by which to filter SV records i.e. 'DEL'; multiple types given space-separated; 'BND' for complex rearrangments with breakends")
   parser.add_argument('-S','--conttype',type=str,nargs='+',help="Complementary type filter, same usage as --T")
   parser.add_argument('--ulim',type=int,help="Provide upper length limit for SVs [bp] (inclusive)")
   parser.add_argument('--llim',type=int,help="Provide lower length limit for SVs [bp] (inclusive)")
   parser.add_argument('--seq',action='store_true',help='Select to write sequence data into single fasta file; requires --out')
   parser.add_argument('--nosize',action='store_true',help='Select to not display the end position and size [bp] of each SV in line')
   parser.add_argument('--circos',action='store_true',help='Circos mode; generates karyotype and link files in current directory for circos for BND-type SVs to show inter-chromosomal events; compatible w/ --C option')
   parser.add_argument('--out',type=str,help='Provide desired name of bed-formatted output file')
   args=parser.parse_args()

   if args.ov:
      if args.c:
         variants=_readfile(args.vcf, cfil=args.contigname, tfil=args.type, vfil=args.conttype)
         _getstats(variants)
      else:
         _contig_overview(args.vcf)
   elif args.circos:
      variants=_readfile(args.vcf, cfil=args.contigname, tfil='BND')
      a=Circos(variants,args.vcf)
   else:
      if args.llim and args.ulim:
         if args.ulim<=args.llim:
            print('Error: Upper limit smaller than lower limit')
            sys.exit()

      variants=_readfile(args.vcf, cfil=args.contigname, tfil=args.type, vfil=args.conttype)
      if args.out:
         [_seq2file(_sortedlist(variants,args.sort), args.out, ulim=args.ulim, llim=args.llim) if args.seq else _list2file(_sortedlist(variants,args.sort), args.out, ulim=args.ulim, llim=args.llim, nosize=args.nosize)]
      else:
         _listputter(_sortedlist(variants,args.sort), ulim=args.ulim, llim=args.llim, nosize=args.nosize)

   print('\n',"--- Run time: {:0.2f} seconds ---".format((time.time() - start_time)))
