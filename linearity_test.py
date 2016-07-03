# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 01:33:41 2016

@author: viherbos
"""

import sys
#from distutils.core import setup
#import pyximport; pyximport.install()
import DBLR as DB
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import FEParam as FEP
from fit_library import *



def linearity_test(base_path, init_point, end_point, step, n_meas, coef, noise,
			  thr_sigma, SPE_DAQ, inc_point, graph_sw):

# Takes linearity DATASET measurements and fits lines to measure the linearity
# of the whole system.
# PARAMETERS:
	#base_path: Base path including the name of the file, number not included (see main)
	#init_point: Width of the first pulse in microseconds
	#end_point:  Width of the last pulse in microseconds
	#step:      Step between pulse in microseconds
	#n_meas:    Number of measurements for every point
	#coef:      Coef for BLR reconstruction
	#noise:     in ADC counts
	#thr_sigma: Number of sigmas of noise for the threshold
	#SPE_DAQ:   Integrated SPE in ADC counts (DAQ_TS = 25 ns)
	#inc_point: First point to start checking saturation. Function will interpolate lines
	#	      for n+1 points showing all the DATASET in order to visually check
	#	      saturation effect.
	#graph_sw:  Graphics switch to show reconstruction signal. Use for debugging purposes only


	########################################################################
	n_points = int((end_point-init_point)/step+1)
	energia_m     = np.zeros(n_points)
	energia_sigma = np.zeros(n_points)
	energia_aux   = np.zeros(n_meas)

	a=0


	for x in range(0,n_points):

		for y in range(0, n_meas):

			path=''.join([base_path,str(x*n_meas+y+1),'.txt'])
			g=pd.read_csv(path)
			f=4096-g.values

			recons,energia_aux[y] = DB.BLR( f.flatten().astype(float), coef,
								   n_sigma = thr_sigma, NOISE_ADC=noise,
								   plot=False)
			if (graph_sw == 1):
				plt.plot(recons)

		energia_m[a]     = np.mean(energia_aux/(SPE_DAQ*FEP.sampling_DAQ()))
		energia_sigma[a] = np.std(energia_aux/(SPE_DAQ*FEP.sampling_DAQ()),ddof=1)
		a=a+1
		if (graph_sw == 1):
			print(energia_m)

	X=np.linspace(init_point,end_point,n_points)

	XI2_r = np.zeros(n_points)

	for l in range(inc_point, n_points+1):

		coeff, perr, XI2_r = line_fit(energia_m[0:l],X[0:l],energia_sigma[0:l],\
								'$\mu$ seconds', 'PE', 'Linear FIT',l,1)

		#Plots fit results for a point by point test
		plt.errorbar(X,energia_m,fmt='b*',yerr=energia_sigma)
		plt.plot(X,line(X,coeff[0],coeff[1]))
		print('CHI_SQUARE_REDUCED =',XI2_r)

	return coeff, perr, XI2_r


def main():
		#PARAMETERS
	init_point = 30
	end_point  = 160
	step       = 10
	n_meas     = 20
	inc_point  = 3
	#PATH
	base_path = 'D:/DATOS_DAC/2052/pmt_0_trace_evt_'

	coeff, perr, XI2_r = linearity_test(base_path, init_point, end_point,
							 step, n_meas, 1.636E-3, 0.75,
						       40, 20.5,inc_point,1)

if __name__ == "__main__":
	main()



