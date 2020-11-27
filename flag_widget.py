import tkinter as tk
import SVread as read
import statistics

def _printstats(self,key,value):
    try:
       low=min(value)
       high=max(value)
       self.longest_SV.append(high)
       avg=sum(value)/len(value)
       med=statistics.median(value)
       #print(f"{key}: {len(value)}  --- MIN: {low}, MAX: {high}, AVG: {avg:.1f}, MED: {med:.0f}")
       self.txt.insert(tk.END, f"{key}: {len(value)}  --- MIN: {low}, MAX: {high}, AVG: {avg:.1f}, MED: {med:.0f}\n")
    except ValueError:
       pass

def _getstats(self):
    self.longest_SV=[]
    vardict=read._labeldict(self.variants)
    [x.sort for x in vardict.values()]
    #print(f'Total number of SVs: {len(self.variants)}')
    _getclen(self)
    self.txt.insert(tk.END, f'\n{len(self.variants)} SVs in total on {self.contig_var.get()} (Length={self.clength} bp);\n')
    [_printstats(self,x,y) for x,y in vardict.items()]
    self.txt.insert(tk.END, '\n')
    
def _getclen(self):
    for x in self.contigs:
        if x.name==self.contig_var.get():
           self.clength=x.length
   

def _contig_varnum(self):
    for line in open(self.vcf,'r'):
        if line[0]!='#':
           for contig in self.contigs:
               if contig.name==line.split()[0]:
                  contig.varnum+=1    