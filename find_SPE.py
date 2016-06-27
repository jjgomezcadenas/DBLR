# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 23:44:40 2016

@author: viherbos
"""

import numpy as np
import pandas as pd

from fit_library import *

#PARAMETERS

#PATH
#base_path = 'D:/DATOS_DAC/spe_1230/2046/pmt_0_trace_evt_'
base_path = '/Users/jjgomezcadenas/Documents/Development/NEXT/data/pmtTraces/pmt_0_trace_evt_'
n_files   = 5000

#INTEGRATING RANGE
start=(43.3*40)
end  =(43.8*40)

#GRAPHICS WINDOW
x_text = 'ADC_COUNTS (LSB)'
y_text = 'Ocurrences'
title_text = 'INTEGRATED SINGLE PHOTOELECTRON VALUE'
bins= 200
n_figure = 1



# Read all files and integrate in the signal region
integral_r = np.zeros(n_files,float)


for x in range(1, n_files):

    path=''.join([base_path,str(x),'.txt'])
    #f = np.loadtxt(path)

    g=pd.read_csv(path)

    f=g.values

    media=np.mean(f)

    f = f[start:end] - media
    integral_r[x-1]=np.sum(f)


a,b = gauss2_fit(integral_r, x_text, y_text, title_text, bins,[0,-20], n_figure,1)


