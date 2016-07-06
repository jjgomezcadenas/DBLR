# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 15:18:42 2016

@author: viherbos
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl

import os


def signal_to_hdf5(in_path, n_points, n_measurements, n_samples, out_path):

	if os.path.exists(out_path):
		os.remove(out_path)

	store = pd.HDFStore(out_path, complevel=9, complib='zlib')
	#store.append('data',pd.DataFrame(np.ones((1,3200), dtype='int32')))
	#store.append('data',pd.DataFrame(np.zeros((1,20), dtype='int32')))
	#f=np.zeros((1,3199),dtype='int32')

	data = pd.DataFrame(dtype='int32',columns=range(0,n_measurements),
				   index=range(0,n_samples))

	for x in range(1, n_points+1):
		# Read all the files in the dataset
		for y in range(1, n_measurements+1):
			path=''.join([in_path,str((x-1)*n_measurements+y),'.txt'])
			# File to read
			g=pd.read_csv(path, names=[str(y-1)], header=None, dtype='int32')
			# DataFrame that stores read data and gives columns names.
			# Each column is a measurement
			data[str(y-1)]=g
			#Store data colum in its place

		store.append("data"+str(x-1),data)
		# Appends new node for every point

	store.close()

def read_hdf5(in_path, point, measurement):

	a = pd.read_hdf(in_path,"data"+str(point))
	b = np.array(a[str(measurement)])

	return b



def main():
	#in_path  = 'D:/DATOS_DAC/spe_1230/2046/pmt_0_trace_evt_'
	in_path  = 'D:/DATOS_DAC/2052/pmt_0_trace_evt_'
	out_path = '2052.h5.z'
	n_points   = 14
	n_measurements = 20
	n_samples = 40000

	signal_to_hdf5(in_path, n_points, n_measurements, n_samples, out_path)

	#res = read_hdf5(out_path, 13, 19)

	#pl.plot(res)
	#pl.show()

if __name__ == "__main__":
	main()