class ToolTip:
    def __init__(self,widget,text,info_var):
        self.text=text
        self.widget=widget
        self.info_var=info_var
        self.widget.bind("<Enter>",self.on_enter)
        self.widget.bind("<Leave>",self.on_leave)
        
    def on_enter(self,event):
        self.info_var.set(self.text)
        
    def on_leave(self,enter):
        self.info_var.set('Hover over widgets for information')