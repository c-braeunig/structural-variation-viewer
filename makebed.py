#!/usr/bin/env python3

import argparse

parser=argparse.ArgumentParser(description='Generate artificial long reads from existing assembly with sliding window approach')
parser.add_argument('--fi',help='Provide index of input fasta file (generate w/ samtools faidx)')
parser.add_argument('--sz',type=int,help='Size of artifical long reads [bp]')
parser.add_argument('--st',type=int,help='Step size')
parser.add_argument('--bed',help='Output file name')
args=parser.parse_args()

def _global(sz,st):
    global size
    global step
    size=sz
    step=st

def _splits(line):
    columns=line.split()
    return [columns[0],int(columns[1])]

def _slide(cont_name,cont_len,outhandle):
    start=0
    stop=size
    while stop<cont_len:
          outhandle.write(f'{cont_name}\t{start}\t{stop}\n')
          start+=step
          stop=start+size
    stop=start+(cont_len-start)
    outhandle.write(f'{cont_name}\t{start}\t{stop}\n')

_global(args.sz,args.st)
with open(args.bed,'w') as outhandle:
     for line in open(args.fi,'r'):
         _slide(*_splits(line),outhandle)
outhandle.close()
