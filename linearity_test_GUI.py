# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 00:31:32 2016

@author: viherbos
"""


from PIL import Image, ImageTk
from Tkinter import Tk, BOTH, W, N, E, S
from Tkinter import Checkbutton, BooleanVar, IntVar, DoubleVar, Spinbox
from Tkinter import StringVar
from ttk import Frame, Button, Label, Entry
import linearity_test as linearity
import matplotlib.pyplot as plt

a=0

class Example(Frame):

	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent
		self.initUI()
		self.centerUI(w=450,h=340)

	def initUI(self):

		self.parent.title("LINEARITY TEST FOR PMT BASES")
		self.pack(fill=BOTH, expand=True)

		self.columnconfigure(0, weight=1)
		#self.rowconfigure(0, weight=1)
		# weight attibute is used to make them growable

#		self.graph_cb   = BooleanVar()
		self.init_point = IntVar()
		self.base_path  = StringVar()
		self.end_point  = IntVar()
		self.step       = IntVar()
		self.n_meas     = IntVar()
		self.inc_point  = IntVar()
		self.coef       = DoubleVar()
		self.noise      = DoubleVar()
		self.thr_sigma  = DoubleVar()
		self.SPE_DAQ    =	 DoubleVar()	
		
		search = Image.open("next_logo.jpg")
		search_temp = search.resize((170, 200), Image.ANTIALIAS)
		search_aux = ImageTk.PhotoImage(search_temp)
		label1 = Label(self, image=search_aux)
		label1.image = search_aux
		label1.grid(row=0, column=0,
				columnspan=10, rowspan=10, sticky=E+W+S+N)

		#Text Box
		self.base_path.set("F:/DATOS_DAC/2052/pmt_0_trace_evt_")
		e1 = Entry(self, textvariable=self.base_path, width=40)
		e1.grid(row=1,column=2, sticky=W, columnspan=5, pady=5)
		e1_label = Label(self, text="DataSet path (including name file)")
		e1_label.grid(row=0,column=2,sticky=W, columnspan=5, pady=5)		
		
		#Spin Boxes
		self.n_meas.set("20")
		sb1 = Spinbox(self, from_=1, to=1000, 
				  width=6, textvariable=self.n_meas)
		sb1.grid(row=3,column=3, sticky=W)
		sb1_label = Label(self, text="Measurements")
		sb1_label.grid(row=2,column=3, padx=0, sticky=W)		
		
		self.step.set("10")
		sb2 = Spinbox(self, from_=10, to=200, 
				  width=6, textvariable=self.step)
		sb2.grid(row=3,column=4, sticky=W)
		sb2_label = Label(self, text="Pulse Width Step")
		sb2_label.grid(row=2,column=4, padx=0, sticky=W)
		
		# INTEGRATION LIMITS
		Integration_label = Label(self, text="INTEGRATION LIMITS",
		                          font = "Verdana 12 bold")
		Integration_label.grid(row=4,column=3, 
						padx=5,
						columnspan = 3, pady=10)
		self.init_point.set("30")
		sb3 = Spinbox(self, from_=1, to=10000, 
				  width=6, textvariable=self.init_point)
		sb3.grid(row=7,column=3, sticky=W)
		sb3_label = Label(self, text="Start (usec)")
		sb3_label.grid(row=6,column=3, padx=0, sticky=W)		
		
		self.end_point.set("160")
		sb4 = Spinbox(self, from_=1, to=10000, 
				  width=6, textvariable=self.end_point)
		sb4.grid(row=7,column=4, sticky=W)
		sb4_label = Label(self, text="End (usec)")
		sb4_label.grid(row=6,column=4, padx=0, sticky=W)

		
		# PARAMETERS
		Integration_label = Label(self, text="PARAMETERS",
		                          font = "Verdana 12 bold")
		Integration_label.grid(row=8,column=3, 
						padx=5,
						columnspan = 3, pady=10)
		self.inc_point.set("3")
		sb5 = Spinbox(self, from_=1, to=100, 
				  width=6, textvariable=self.inc_point)
		sb5.grid(row=11,column=3, sticky=W)
		sb5_label = Label(self, text="First point")
		sb5_label.grid(row=10,column=3, padx=0, sticky=W)
		
		self.coef.set("1.636E-3")
		e6 = Entry(self, width=10, textvariable=self.coef)
		e6.grid(row=11,column=4, sticky=W)
		e6_label = Label(self, text="DBLR Coef")
		e6_label.grid(row=10,column=4, sticky=W)
				
		self.noise.set("0.75")
		e7 = Entry(self, width=10, textvariable=self.noise)
		e7.grid(row=11,column=5, sticky=W)
		e7_label = Label(self, text="Noise (LSB)")
		e7_label.grid(row=10,column=5, sticky=W)

		self.thr_sigma.set("40")
		e8 = Entry(self, width=10, textvariable=self.thr_sigma)
		e8.grid(row=13,column=3, sticky=W)
		e8_label = Label(self, text="Threshold")
		e8_label.grid(row=12,column=3, sticky=W)
		
		self.SPE_DAQ.set("20.5")
		e9 = Entry(self, width=10, textvariable=self.SPE_DAQ)
		e9.grid(row=13,column=4, sticky=W)
		e9_label = Label(self, text="SPE (LSB)")
		e9_label.grid(row=12,column=4, sticky=W)		
		
		
#		#Check buttons
#		cb1 = Checkbutton(self, text="MultiGraph Output", variable=self.graph_cb					)
#		cb1.select()
#		cb1.grid(row=7,column=6, sticky=W)

		
		# Main buttons
		obtn = Button(self, text="GO!!", command=self.linearity_f)
		obtn.grid(row=14, column=4, sticky=E, pady=10)

		cbtn = Button(self, text="Quit", command=self.quit)
		cbtn.grid(row=14, column=5, sticky=E, pady=10)

		hbtn = Button(self, text="Help")
		hbtn.grid(row=14, column=0, sticky=W, pady=10)

	def linearity_f(self):
		
#		global a
#		
#		if (self.graph_cb.get()==True): 
#			b=a
#			a=a+1
#		else: 
#			b=0				
		
		linearity.linearity_test(base_path=self.base_path.get(),
				              init_point=self.init_point.get(), 
				              end_point=self.end_point.get(), 
				              step=self.step.get(), 
				              n_meas=self.n_meas.get(), 
				              inc_point=self.inc_point.get(), 
				              coef=self.coef.get(),
				              noise=self.noise.get(),
					        thr_sigma=self.thr_sigma.get(),
					        SPE_DAQ=self.SPE_DAQ.get(),
					        graph_sw=0
						)
		
		
	def centerUI(self,w,h):
		sw = self.parent.winfo_screenwidth()
		sh = self.parent.winfo_screenheight()
		x  = (sw-w)/2
		y  = (sh-h)/2
		self.parent.geometry('%dx%d+%d+%d' % (w,h,x,y))

def main():

	root = Tk()
	root.resizable(width=False, height=False)
	app = Example(root)
	root.mainloop()


if __name__ == '__main__':
    main()