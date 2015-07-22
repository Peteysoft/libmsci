#!/usr/bin/python

import sys
import Tkinter as tk
import Image

#imagefile="/media/Acer/cygwin/home/Petey/homepages/libagf_webpage/pics/N_Netherlands_waterways800.gif"
#imagefile="image2.pgm"
#imagefile="/media/Acer/cygwin/home/Petey/data/landsat/LE70160282011169EDC00_B1.TIF"

def mouseclick(event):
    canvas=event.widget
    x=canvas.canvasx(event.x)
    y=canvas.canvasy(event.y)
    print x, y, event.num-1
    fs.write(str(x)+" "+str(y)+" "+str(event.num-1)+"\n")

class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.grid(sticky=tk.N+tk.S+tk.E+tk.W)
        self.createWidgets()

    def mousewheel(self, event):
        self.screen.yview_scroll(-1*(event.delta/120), "units")

    def createWidgets(self):
        top=self.winfo_toplevel()
        top.rowconfigure(0, weight=1)
        top.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        #self.im=Image.open(imagefile)
        #print (self.im.size[0]*2, self.im.size[1]*2)
        #self.im.resize((self.im.size[0]*2, self.im.size[1]*2))
        #self.image=tk.PhotoImage(self.im)
        self.image=tk.PhotoImage(file=imagefile)
        self.screen = tk.Canvas(self, width=800, height=600, scrollregion=(0, 0, self.image.width(), self.image.height()), confine=1, cursor="crosshair")
        self.screen.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        self.hbar=tk.Scrollbar(self, orient=tk.HORIZONTAL, command=self.screen.xview)
        self.hbar.grid(row=1, column=0, sticky=tk.E+tk.W)
        self.vbar=tk.Scrollbar(self, orient=tk.VERTICAL, command=self.screen.yview)
        self.vbar.grid(row=0, column=1, sticky=tk.N+tk.S)
        self.screen.config(xscrollcommand=self.hbar.set, yscrollcommand=self.vbar.set)
        self.screen.create_image(0,0, anchor=tk.NW, image=self.image)
        self.screen.bind_all('<MouseWheel>', self.mousewheel)
        self.screen.bind('<Button-1>', mouseclick)
        self.screen.bind('<Button-2>', mouseclick)
        self.screen.bind('<Button-3>', mouseclick)
        #self.screen.bind('<space>', mouseclick)
        #self.screen.bind('<Enter>', mouseclick)
        #self.screen.bind('<Tab>', mouseclick)

imagefile=sys.argv[1]
fs=open(sys.argv[2], "w")        
app = Application()
app.master.title('Sample application')
app.mainloop()

fs.close()

