# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 01:33:41 2016

@author: viherbos
"""
from PlotUtil import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from MauDeconv2 import MauDeconv
from fit_library import *
import DBLR as DB

#PARAMETERS
init_point = 30
end_point  = 30
step       = 10
n_meas = 1
n_points = (end_point-init_point)/step+1

#PATH
#base_path = 'F:/DATOS_DAC/2052/pmt_0_trace_evt_'
base_path = '/Users/jjgomezcadenas/Documents/Development/NEXT/data/pmtTraces/pmt_0_trace_evt_'

###########################################################################
energia_m     = np.zeros(n_points)
energia_sigma = np.zeros(n_points)
energia_aux   = np.zeros(n_meas)

a=0


for x in range(0,n_points):

	for y in range(0, n_meas):

		path=''.join([base_path,str(x*n_meas+y+1),'.txt'])
		g=pd.read_csv(path)
		f=g.values
		fs = np.array(4096-f)
		fs2 = np.reshape(fs,fs.shape[0])
		# recons, energia_aux[y] = MauDeconv([(1)], 4096-f, 1.636E-3, \
		# 				n_sigma=150, MAU_WindowSize=256, SPE=20.5)

		signal_daq = np.zeros(len(f), dtype=np.double)

		print "shape of f =", f.shape, "shape of fs2 =", fs2.shape
		
		signal_t = np.arange(0.0, len(signal_daq)*1., 1., dtype=np.double)
		signal_daq[:] = fs2

		signal_rec = DB.BLR(signal_t, signal_daq, 1.636E-3, 
									    n_sigma = 2, NOISE_ADC=0.7)

		#energia_aux[y] = rec

	#plt.plot(recons)
	plot_signal(signal_t,signal_rec, 
                title = 'signal', signal_start=20000, signal_end=21000, units='ns')
	plt.show()	
	
	energia_m[a]     = energia_aux.mean()
	energia_sigma[a] = energia_aux.std()
	a=a+1 
	print energia_m

#X=np.linspace(init_point,end_point,n_points)
#
#aXI2_r = np.zeros(n_points)
#
#for l in range(3, 5):
#		
#	coeff, perr, XI2_r = line_fit(energia_m[0:l],X[0:l],energia_sigma[0:l],\
#							'uSeconds', 'PE', 'Linear FIT',l,0)
#	
#	#Plots fit results for a point by point test
#	plt.errorbar(X,energia_m,fmt='b*',yerr=energia_sigma)
#	plt.plot(X,line(X,coeff[0],coeff[1]))
#	print 'CHI_SQUARE_REDUCED =',XI2_r







