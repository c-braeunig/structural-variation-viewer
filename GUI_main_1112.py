import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg

import re
import SVread as read
import SVsee as see
import SVpipe as pipe
import flag_widget as fw
import tool_tip as tp



class Contig:
      def __init__(self,name,length):
          self.name=name.rstrip().lstrip()
          self.length=length
          self.varnum=0
      def add2list(self,targetlist):
          targetlist.append(self)

class CustomToolbar(NavigationToolbar2TkAgg):
    def __init__(self,canvas_,parent_):
        self.toolitems = (
            ('Home', 'Reset view', 'home', 'home'),
            (None, None, None, None),
            ('Pan', 'Pan', 'move', 'pan'),
            ('Zoom', 'Zoom', 'zoom_to_rect', 'zoom'),
            (None, None, None, None),
            ('Save', 'Save figure', 'filesave', 'save_figure')
            )
        NavigationToolbar2TkAgg.__init__(self,canvas_,parent_)

class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self,master)

        self.initUI()
        master.wm_title("Visualize")
        master.geometry("975x400")
        #self.master.tk.call('wm', 'iconphoto', 
        #                    self.master._w,
        #                    tk.PhotoImage(file='leopardV.png'))

    def initUI(self):
        self.main_menu()
        
        self.subframe_right=tk.Frame(self.master)
        self.subframe_right.grid(row=0,column=1,padx=5,pady=3,sticky='nsew')
        
        self.p = ttk.Panedwindow(self.subframe_right, orient=tk.HORIZONTAL)
        self.q = ttk.Panedwindow(self.subframe_right, orient=tk.VERTICAL)

        self.draw_pane=ttk.Labelframe(self.p,text='Visualize',
                                      width=75,height=75,labelanchor='n')
        self.contig_pane = ttk.Labelframe(self.p, text='Select contig',
                                          width=150, height=75)
        self.p.add(self.contig_pane,weight=1)
        self.p.add(self.draw_pane,weight=0)
        self.p.grid(row=0,column=1,sticky='new')
        
        self.parameter_pane = ttk.Labelframe(self.q, text='Main parameters',
                                             width=150, height=200)
        self.annotation_pane = ttk.Labelframe(self.q, text='Annotation parameters',
                                              width=150, height=80)
        self.info_pane=ttk.Labelframe(self.q, text='Info',
                                      width=150, height=150)
        self.q.add(self.parameter_pane)
        self.q.add(self.annotation_pane)
        self.q.add(self.info_pane)
        self.q.grid(row=1,column=1,sticky='new')
        self.info_pane.columnconfigure(0,weight=5)
         
        self.info_box()
        self.text_frame()        
        self.contig_menu()
        self.type_menu()
        self.draw_button()
        self.heat_menu()
        self.limit_menu()
        self.region_menu()
        self.anno_check()
        self.anno_type_menu()
        self.anno_limit_menu()
        
    def text_frame(self):
        self.subframe_left = ttk.Frame(self.master,width=500,height=300)
        self.txt = tk.Text(self.subframe_left)
        
        self.subframe_left.grid(row=0,column=0,sticky='nw')
        self.txt.grid(row=0,column=0)

    def main_menu(self):
        self.annotation_file=None
        
        menubar = tk.Menu(self)
        
        fileMenu=tk.Menu(self,tearoff=0)
        fileMenu.add_command(label="Variant file", command=self.onOpen)
        fileMenu.add_command(label='Annotation file', command=self.open_annotation)
        fileMenu.add_command(label='Unload annotation file',command=self.unload)
        menubar.add_cascade(label="Load", menu=fileMenu)
        
        menubar.add_cascade(label='Close',command=self.quit_all)
        
        self.master.config(menu=menubar)
        
    def quit_all(self):
        try:
           self.plot.close('all')
           self.quit()
        except AttributeError:
           self.quit()
        
    def onOpen(self):
        ftypes = [('Variant Call Format','.vcf'),('All files', '*')]
        dlg = filedialog.Open(self, filetypes = ftypes, initialdir=".")
        fl = dlg.show()
        if fl != '':
            self.vcf=fl
            self.readvcf()
            self.txt.insert(tk.END, f'{fl} loaded...\n')
            self.txt.insert(tk.END, f'{len(self.contigs)} contigs found...\n')
            self.update_option_menu()
            
    def info_box(self):
        self.info_var=tk.StringVar(self.info_pane)
        self.info_var.set('Hover over widgets for information')
        self.info_box=tk.Label(self.info_pane,
                               textvariable=self.info_var,
                               height=5,width=25,
                               justify=tk.CENTER)
        self.info_box.grid(row=0,column=0,sticky='nsew',pady=3,padx=5)
        
        self.info_box_var=self.info_var
    
    def contig_menu(self):
        self.contig_var=tk.StringVar(self.contig_pane)
        self.default=Contig('------',0)
        self.contigs=[self.default]
        
        self.contig_var.set(self.default.name)
        
        self.over_drop=tk.OptionMenu(self.contig_pane,
                                     self.contig_var,
                                     *[x.name for x in self.contigs])       
        self.over_drop.config(width=30)
        self.over_drop.grid(row=0,column=0,padx=3,pady=3,sticky='nw')
        self.contig_var.trace('w',self.flagstat)
        
        tp.ToolTip(widget=self.over_drop,
                   text='Select contig to visualize;\n Only one contig can be shown at a time',
                   info_var=self.info_box_var)
        
    def update_option_menu(self):
        menu = self.over_drop["menu"]
        menu.delete(0, "end")
        for string in self.contigs:
            menu.add_command(label=f'{string.name} ({string.varnum})',
                             command=lambda value=string.name: self.contig_var.set(value))
    
    def flagstat(self,*args):
        if self.contig_var.get()==self.default.name:
           self.info_var.set('Please first load a vcf file\n before trying to select a contig\n')
        else:
           self.variants=read._readfile(self.vcf,cfil=[self.contig_var.get()])
           fw._getstats(self)
           self.update_type_menu()
           self.set_default_length()
           self.set_default_limit()
           try:
              self.annotation_overview()
              self.update_anno_type_menu()
           except AttributeError:
              pass
           except TypeError:
              pass

    def heat_menu(self):
        self.heat_var=tk.IntVar(self.parameter_pane)
        self.heat_entry=tk.Entry(self.parameter_pane,state=tk.DISABLED,
                                 textvariable=self.heat_var,width=7)
        self.heat_var.set(0)

        self.heat_check=tk.BooleanVar(self.parameter_pane)
        self.heat_check.set(False)
        self.heat_checkButton=tk.Checkbutton(self.parameter_pane,text='Heat map',
                                             var=self.heat_check,
                                             command=self.enable_heat_entry)
                
        self.heat_checkButton.grid(row=3,column=0,padx=3,pady=0,sticky='w')
        self.heat_entry.grid(row=3,column=1)
        
        tp.ToolTip(widget=self.heat_checkButton,
                   text='Check to enable heatmap for the density of SVs.\n Enter distance (in bp) over which density is calculated',
                   info_var=self.info_box_var)

    def enable_heat_entry(self):
        if self.heat_check.get():
           self.heat_entry.configure(state=tk.NORMAL)
           self.heat_var.set(0)
        else:
           self.heat_entry.delete(0,'end')
           self.heat_entry.configure(state=tk.DISABLED)
           self.heat_var.set(0)

    def type_menu(self):
        self.typeVerbose=tk.BooleanVar(self.parameter_pane)
        self.typeVerbose.set(False)
        
        self.typeButton=tk.Menubutton(self.parameter_pane,text='Type selection',
                                      relief=tk.RAISED)
        self.typeButton.grid(row=0,column=0,sticky='w',columnspan=2)
        self.VerboseButton=tk.Checkbutton(self.parameter_pane,text='Verbose selection',
                                          variable=self.typeVerbose)
        self.VerboseButton.grid(row=0,column=1)
        
        tp.ToolTip(widget=self.typeButton,
                   text='Select SV type(s) to visualize.\n No selection means all types will be shown.\n Checking "Verbose selection" means all types\n but the one(s) selected will be shown',
                   info_var=self.info_box_var)
                
    def update_type_menu(self):
        self.types=list(set([x.type for x in self.variants]))
        self.typeButton.menu  =  tk.Menu ( self.typeButton, tearoff = 0 )
        self.typeButton["menu"]  =  self.typeButton.menu

        self.type_selection=[]
        self.typeButton.menu.delete(0, "end")
        for x in self.types:
            exec(f'self.{x}_var=tk.BooleanVar(self.master)')
            exec(f'self.{x}_var.set(False)')
            exec(f'self.typeButton.menu.add_checkbutton(label="{x}",variable=self.{x}_var)')

    def limit_menu(self):
        self.lower=tk.IntVar(self.parameter_pane)
        self.upper=tk.IntVar(self.parameter_pane)
        
        self.lower_entry=tk.Entry(self.parameter_pane,state=tk.DISABLED,
                                  textvariable=self.lower,width=10)
        self.upper_entry=tk.Entry(self.parameter_pane,state=tk.DISABLED,
                                  textvariable=self.upper,width=10)
        
        self.limit_check=tk.BooleanVar(self.parameter_pane)
        self.limit_check.set(False)
        self.limit_checkButton=tk.Checkbutton(self.parameter_pane,text='Element size limits',
             var=self.limit_check,
             command=self.enable_limit_entry)
        
        self.limit_checkButton.grid(row=1,column=0,sticky='w')
        self.lower_entry.grid(row=1,column=1)
        self.upper_entry.grid(row=1,column=2)
        
        tp.ToolTip(widget=self.limit_checkButton,
                   text='Select bp length limits \n for SVs to be shown',
                   info_var=self.info_box_var)
                       
    def enable_limit_entry(self):
        if self.limit_check.get():
           self.lower_entry.configure(state=tk.NORMAL)
           self.upper_entry.configure(state=tk.NORMAL)
        else:
           self.lower_entry.configure(state=tk.DISABLED)
           self.upper_entry.configure(state=tk.DISABLED)
           self.set_default_limit()
           
    def set_default_limit(self):
        #if self.lower.get()==0 and self.upper.get ##try to not keep resetting if 
        ##limits have already been selected...same for anno limits and regions
        if not self.limit_check.get():
           self.lower.set(0)
           try:
              self.upper.set(max(self.longest_SV))
           except ValueError:
              self.upper.set(0)
           except AttributeError:
              self.info_var.set('Please first load a vcf file\n before trying to set limits\n')
           
    def anno_limit_menu(self):
        self.a_lower=tk.IntVar(self.annotation_pane)
        self.a_upper=tk.IntVar(self.annotation_pane)
    
        self.a_lower_entry=tk.Entry(self.annotation_pane,state=tk.DISABLED,
                                    textvariable=self.a_lower,width=10)
        self.a_upper_entry=tk.Entry(self.annotation_pane,state=tk.DISABLED,
                                    textvariable=self.a_upper,width=10)
        
        self.a_limit_check=tk.BooleanVar(self.annotation_pane)
        self.a_limit_check.set(False)
        self.a_limit_checkButton=tk.Checkbutton(self.annotation_pane,
                                                text='Element size limits',
                                                var=self.a_limit_check,
                                                command=self.enable_a_limit_entry)
        
        self.a_limit_checkButton.grid(row=2,column=0,sticky='w')
        self.a_lower_entry.grid(row=2,column=1)
        self.a_upper_entry.grid(row=2,column=2)
        self.annotation_pane.columnconfigure(0,weight=0)
        self.annotation_pane.columnconfigure(1,weight=1)
        
        tp.ToolTip(widget=self.a_limit_checkButton,
                   text='Select bp length limits for annotation \n features to be shown',
                   info_var=self.info_box_var)
        
    def enable_a_limit_entry(self):
        if self.a_limit_check.get():
           self.a_lower_entry.configure(state=tk.NORMAL)
           self.a_upper_entry.configure(state=tk.NORMAL)
        else:
           self.a_lower_entry.configure(state=tk.DISABLED)
           self.a_upper_entry.configure(state=tk.DISABLED)
           self.set_default_a_limit()
           
    def set_default_a_limit(self):
        if not self.a_limit_check.get():
           self.a_lower.set(0)
           try:
              self.a_upper.set(self.longest_anno)
           except ValueError:
              self.a_upper.set(0)
           except AttributeError:
              self.info_var.set('Please first load an annotation file\n before trying to set limits\n')
           
    def anno_check(self):
        self.anno_check=tk.BooleanVar(self.annotation_pane)
        self.anno_check.set(False)
        self.anno_checkButton=tk.Checkbutton(self.annotation_pane,text='Draw annotations',
                                             var=self.anno_check)            
        self.anno_checkButton.grid(row=1,column=0,padx=3,pady=0,sticky='w')
        
        tp.ToolTip(widget=self.anno_checkButton,
                   text='If annotation file is loaded,\n annotations are drawn by default.\n To disable, toggle checkbox',
                   info_var=self.info_box_var)
    
    def region_menu(self):
        self.lreg=tk.IntVar(self.parameter_pane)
        self.ureg=tk.IntVar(self.parameter_pane)
        
        self.lreg_entry=tk.Entry(self.parameter_pane,state=tk.DISABLED,
                                 textvariable=self.lreg,width=10)
        self.ureg_entry=tk.Entry(self.parameter_pane,state=tk.DISABLED,
                                 textvariable=self.ureg,width=10)
        
        self.reg_check=tk.BooleanVar(self.parameter_pane)
        self.reg_check.set(False)
        self.reg_checkButton=tk.Checkbutton(self.parameter_pane,text='Region selection',
             var=self.reg_check,
             command=self.enable_reg_entry)
        
        self.reg_checkButton.grid(row=2,column=0,sticky='w')
        self.lreg_entry.grid(row=2,column=1)
        self.ureg_entry.grid(row=2,column=2)
        
        tp.ToolTip(widget=self.reg_checkButton,
                   text='Select bp coordinates between which too show SVs',
                   info_var=self.info_box_var)
                
    def set_default_length(self):
        if not self.reg_check.get():
           for x in self.contigs:
               if x.name==self.contig_var.get():
                  self.clen=x.length
                  self.lreg.set(0)
                  self.ureg.set(x.length)
    
    def enable_reg_entry(self):
        if self.reg_check.get():
           self.lreg_entry.configure(state=tk.NORMAL)
           self.lreg.set(0)
           self.ureg_entry.configure(state=tk.NORMAL)
        else:
           self.lreg_entry.configure(state=tk.DISABLED)
           self.ureg_entry.configure(state=tk.DISABLED)        
          
    def open_annotation(self):
        ftypes = [('GFF files','.gff'),('BED files', '.bed'),('BLAST output','*')]
        dlg = filedialog.Open(self, filetypes = ftypes)
        fl = dlg.show()
        if fl!='':
           self.annotation_file=fl
           self.txt.insert(tk.END, f'{fl} loaded for annotation\n\n')
           self.anno_check.set(True)
           if self.contig_var.get()!=self.default.name:
              self.annotation_overview()
    
    def unload(self):
        try:
           if self.annotation_file:
              self.txt.insert(tk.END,f'{self.annotation_file} unloaded\n\n')
           del(self.annotation_file)
           self.anno_check.set(False)
           self.a_limit_check.set(False)
           self.a_lower.set(0)
           self.a_upper.set(0)
           self.a_typeVerbose.set(False)
           self.a_typeButton.menu.delete(0, "end")
           del(self.annolist)
        except AttributeError:
           self.info_box_var.set('No file to unload')
           
    def annotation_overview(self):
        contig=self.contig_var.get()
        self.ext=re.split('\.',self.annotation_file)[-1]
        self.annolist=see._tupperware(self.annotation_file,contig,type=self.ext)
        if self.annolist:
           self.longest_anno=max([x.length for x in self.annolist])
           self.set_default_a_limit()
           self.update_anno_type_menu()
           self.txt.insert(tk.END,f'Found {len(self.annolist)} annotations on {contig}\n\n')
        else:
           self.txt.insert(tk.END,'No annotations found; make sure vcf and annotation files match\n\n')
           
    def anno_type_menu(self):
        self.a_typeVerbose=tk.BooleanVar(self.annotation_pane)
        self.a_typeVerbose.set(False)
        
        self.a_typeButton=tk.Menubutton(self.annotation_pane,text='Type selection',
                                      relief=tk.RAISED)
        self.a_typeButton.grid(row=0,column=0,sticky='w',columnspan=2)
        self.a_VerboseButton=tk.Checkbutton(self.annotation_pane,text='Verbose selection',
                                          variable=self.a_typeVerbose)
        self.a_VerboseButton.grid(row=0,column=1)
        
        tp.ToolTip(widget=self.a_typeButton,
                   text='Select annotation type(s) to visualize.\n No selection means all types will be shown.\n Checking "Verbose selection" means all types\n but the one(s) selected will be shown\n (Only works with GFF files)',
                   info_var=self.info_box_var)
                
    def update_anno_type_menu(self):
        if self.ext=='gff':
           self.a_types=list(set([x.type for x in self.annolist]))
           self.a_typeButton.menu  =  tk.Menu ( self.a_typeButton, tearoff = 0 )
           self.a_typeButton["menu"]  =  self.a_typeButton.menu

           self.a_type_selection=[]
           self.a_typeButton.menu.delete(0, "end")
           for x in self.a_types:
               exec(f'self.{x}_var=tk.BooleanVar(self.master)')
               exec(f'self.{x}_var.set(False)')
               exec(f'self.a_typeButton.menu.add_checkbutton(label="{x}",variable=self.{x}_var)')
        else:
           self.a_types=None

    def readvcf(self):
        self.contigs=[]
        cont_name=re.compile('<ID=.*,')
        cont_len=re.compile(',length=[0-9]*>')
        with open(self.vcf,'r') as self.inhandle:
             for line in self.inhandle:
                 if line[:8]=='##contig':
                    name=cont_name.search(line)[0].split('=')[1].rstrip(',')
                    length=int(cont_len.search(line)[0].split('=')[1].rstrip('>'))
                    a=Contig(name,length)
                    a.add2list(self.contigs)
        fw._contig_varnum(self)
        self.contigs.sort(key=lambda x: x.length, reverse=True)

    def draw_button(self):
        self.drawButton=tk.Button(self.draw_pane,text='Draw',command=self.draw_window)
        self.drawButton.grid(row=0,column=0)
        tp.ToolTip(widget=self.drawButton,
                   text='Click to visualize selected contig\n and corresponding SVs based on the set parameters',
                   info_var=self.info_box_var)
        
    def filter_variants(self,input,verbose_variable,**kwargs):
        input_new=[]
        lower_limit=kwargs.get('lower',None)
        upper_limit=kwargs.get('upper',None)
        tfil=kwargs.get('tfil',None)
        def filter_variants_by_type(input,tfil,verbose_variable):
            if verbose_variable.get():
               input_new=[x for x in input if x.type not in tfil]
               return input_new
            else:
               input_new=[x for x in input if x.type in tfil]
               return input_new
        try:
           for x in input:
               if lower_limit<=x.length<=upper_limit:
                  input_new.append(x)
        except AttributeError:
           try:
              for x in input:
                  if lower_limit<=x.length:
                     input_new.append(x)
           except AttributeError:
              try:
                 for x in input:
                     if x.length<=upper_limit:
                        input_new.append(x)
              except AttributeError:
                 input_new=input[:]
        if tfil:
           input_new=filter_variants_by_type(input_new,tfil,verbose_variable)
        return input_new
      
    def retrieve_filters(self):
        self.cfil=self.contig_var.get()
        self.tfil=[]
        for x in self.types:
            eval_statement=f'self.{x}_var.get()'
            if eval(eval_statement):
               self.tfil.append(x)
        self.heat=self.heat_var.get()
        self.anno=self.anno_check.get()
        self.llim=self.lower.get()
        self.ulim=self.upper.get()
        self.reg=[self.lreg.get(),self.ureg.get()]
        self.a_llim=self.a_lower.get()
        self.a_ulim=self.a_upper.get()
        
    def pop_up_save_menu(self):
        self.pop_up_save=tk.Menu(self.window, tearoff=0)
        self.pop_up_save.add_command(label='Save figure',command=self.save_figure)
        
    def pop_up(self, event):
        try:
            self.pop_up_save.tk_popup(event.x_root, event.y_root, 0)
        finally:
            self.pop_up_save.grab_release()
            
        
    def save_annotations(self):
        ftypes = [('Text file','.txt'),
                  ('All files', '*')]
        default_file=f'Annotations_{self.contig_var.get()}'
        name=filedialog.asksaveasfile(initialfile=default_file, initialdir="/",
                                      defaultextension='.txt', title="Save file",
                                      filetypes = ftypes)
        if name!='':
           with open(name.name,'w') as output:
                output.write(self.anno_text.get(1.0,tk.END))
        
    def onpick(self,event):
        def get_anno_window(self):       
            self.anno_window=tk.Toplevel(self.window)     
            self.anno_text=tk.Text(self.anno_window)
            
            save_button=tk.Button(self.anno_window,text='Save',command=self.save_annotations)
            save_button.grid(row=0,column=0,sticky='w')
            
            self.anno_text.grid(row=1,column=0,sticky='nsew')
            self.anno_text.insert(tk.END,f'#Variant file: {self.vcf}\n')
            self.anno_text.insert(tk.END,f'#Annotation file: {self.annotation_file}\n')
            self.anno_text.insert(tk.END,f'#Chromosome/Contig: {self.contig_var.get()}\n')
        def get_anno(self):
            if not self.anno_window.winfo_exists():
               get_anno_window(self) 
               
            for x in self.new_anno:
                if x.start==self.xdata and x.length==self.width and x.info==self.label:
                   self.anno_text.insert(tk.END,f'\n#Start: {x.start}, Length: {x.length}\n{x.info}')

        thispatch = event.artist
        self.xdata=thispatch.get_x()
        self.width=thispatch.get_width()
        self.label=thispatch.get_label()
        try:
           get_anno(self)
        except AttributeError:
           get_anno_window(self)
           get_anno(self)
        
        
    def draw_window(self):
        _continue=False
    
        def annotation_types(self):
            self.a_tfil=[]
            if self.a_types:
               for x in self.a_types:
                   eval_statement=f'self.{x}_var.get()'
                   if eval(eval_statement):
                      self.a_tfil.append(x)
        try:
           self.retrieve_filters()
           self.new_variants=self.filter_variants(self.variants,
                                                  self.typeVerbose,
                                                  lower=self.llim,
                                                  upper=self.ulim,
                                                  tfil=self.tfil)
           self.txt.insert(tk.END,f'{len(self.new_variants)} SVs to be drawn\n') 
           try:
              annotation_types(self)
              self.new_anno=self.filter_variants(self.annolist,
                                                 self.a_typeVerbose,
                                                 lower=self.a_llim,
                                                 upper=self.a_ulim,
                                                 tfil=self.a_tfil)
              if self.new_anno:
                 self.txt.insert(tk.END,f'{len(self.new_anno)} annotations to be drawn\n')
              else:
                 self.txt.insert(tk.END,'0 annotations to be drawn; check annotation file and try again\n')
              _continue=True
           except ValueError:
              pass
              _continue=True
        except ValueError:
           self.info_var.set('Please first load a vcf file and select a contig\n')
        
        if _continue:
           self.window = tk.Toplevel()
           #self.window.tk.call('wm', 'iconphoto', self.window._w, 
           #                    tk.PhotoImage(file='leopardV.png'))
        
           self.draw_plot()
           self.window_title=self.plot.ax1.get_title()
           self.window.title(self.window_title)
        
           try:   
              self.canvas = FigureCanvasTkAgg(self.plot.fig, self.window)
              self.canvas.draw()
              self.toolbar = CustomToolbar(self.canvas,self.window)
              self.toolbar.update()
              self.toolbar.pack(side=tk.TOP, fill=tk.X)
              self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH)
              self.canvas.mpl_connect('pick_event', self.onpick)
           except UnboundLocalError:
              print('Unable to draw data; review settings and try again\n')
              self.quit()
           
    def draw_plot(self):
        #if self.annotation_file and self.heat:
        if self.anno and self.heat:
           self.plot = see.Plot(self.cfil,self.clen,reg=self.reg, anno=self.annolist, heat=self.heat)
           pipe._vdraw2(self.new_variants,self.plot)
           pipe._annodraw2(self.new_anno,self.ext,self.plot)
           see._heatmap(self.new_variants,self.clen,self.heat,self.plot,reg=self.reg)
        #elif self.annotation_file:
        elif self.anno:
           self.plot = see.Plot(self.cfil,self.clen,reg=self.reg, anno=self.annolist, heat=self.heat)
           pipe._vdraw2(self.new_variants,self.plot)
           pipe._annodraw2(self.new_anno,self.ext,self.plot)
        elif self.heat: 
           self.plot = see.Plot(self.cfil,self.clen,reg=self.reg, heat=self.heat)
           pipe._vdraw2(self.new_variants,self.plot)
           see._heatmap(self.new_variants,self.clen,self.heat,self.plot,reg=self.reg)
        else:
           self.plot = see.Plot(self.cfil,self.clen,reg=self.reg)
           pipe._vdraw2(self.new_variants,self.plot)
        

    
if __name__=='__main__':
   root = tk.Tk()
   root.resizable(False, False) 
   ex = Application(root)
   root.mainloop()