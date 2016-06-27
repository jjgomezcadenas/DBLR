#from Util import *
#from PlotUtil import *
#import FEE
#import SPE as SP
#import FEParam as FP
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal as SGN

plot_MAU = False
plot_pulse_f = True
plot_tail_f = False
plot_trigger_f = False
plot_signal_i = False
plot_acum = False
plot_signal_daq = True
plot_signal_r = True

def MauDeconv(signal_t_daq, signal_daq, coef, thr = 150,\
		   MAU_WindowSize = 512, SPE = 20.5):
	"""
	Deconvolution offline of the DAQ signal using a MAU
	moving window-average filter of a vector data
	y(n) = (1/WindowSize)(x(n) + x(n-1) + ... + x(n-windowSize))
	in a filter operation filter(b,a,x):
	b = (1/WindowSize)*ones(WindowSize) = (1/WS)*[1,1,1,...]: numerator
	a = 1 : denominator
	y = filter(b,a,x)
	y[0] = b[0]*x[0] = (1/WS) * x[0]
	y[1] = (1/WS) * (x[0] + x[1])
	y[WS-1] = mean(x[0:WS])
	y[WS] = mean(x[1:WS+1])
	and so on
	"""

	MAU = np.zeros(len(signal_daq))
	acum = np.zeros(len(signal_daq))
	signal_r = np.zeros(len(signal_daq))

	pulse_f = np.zeros(len(signal_daq))
	tail_f = np.zeros(len(signal_daq))
	trigger_f = np.zeros(len(signal_daq))

	signal_i = np.copy(signal_daq) #uses to update MAU while procesing signal


	thr2 = 50
	thr_cmp = 2*thr  # to follow baseline tail
	thr_tr = 5. # to conclude BLR when signal_deconv = signal_raw

	#MAU_WindowSize = 40 # provisional
	nm = MAU_WindowSize

 	B_MAU       =   np.ones(nm)

	# print """
	# Mau Window Size = %d
	# number of sigmas = %d
	# noise (adc counts) =%7.2f
	# Mau threshold = %7.2g
	# """%(nm, n_sigma, FP.NOISE_ADC, thr)

# 	print B_MAU

 	#wait()

# 	MAU averages the signal in the initial tranch
# 	 allows to compute the baseline of the signal

    	aux = SGN.lfilter(B_MAU,1, signal_daq[0:nm])
	MAU[0:nm] = aux.transpose()

	# plot_signal2(signal_t_daq[0:nm]/ns,MAU[0:nm],
 #                title = 'MAU: init',
 #                signal_start=0, signal_end=nm*FP.time_bin*ns,
 #                signal_min=-0.5, signal_max=0.5,
 #                units='adc counts')

	# print """
	# MAU: value for n = nm = %7.2f rms = %7.2f
	# """%(MAU[nm], mau_rms)


	acum[nm-1] =  MAU[nm-1]
	BASELINE = MAU[nm-1]


#----------

# While MAU inits BLR is switched off, thus signal_r = signal_daq

	signal_r[0:nm] = signal_daq[0:nm].transpose()

	pulse_on=0
	wait_over=0

	offset = 0

	# MAU has computed the offset using nm samples
	# now loop until the end of DAQ window
	for k in xrange(nm,len(signal_daq)):

		trigger_line = MAU[k-1] + thr

		# follower variables for monitor purposes
		pulse_f[k] = pulse_on + 0.5
		tail_f[k] = wait_over + 0.5
		trigger_f[k] = trigger_line

		# condition: raw signal raises above trigger line and we are not in the tail
		# (wait_over == 0)
		if signal_daq[k] > trigger_line and wait_over == 0:

			# if the pulse just started pulse_on = 0.
			# In this case compute the offset as value of the MAU before pulse starts (at k-1)

			if pulse_on == 0: # pulse just started
				offset = MAU[k-1]  #offset computed as the value of MAU before pulse starts
				pulse_on = 1

			#Freeze the MAU
			MAU[k] = MAU[k-1]

			signal_i[k] =MAU[k-1]  #signal_i follows the MAU

			#update recovered signal, correcting by offset
			acum[k] = acum[k-1] + signal_daq[k] - offset
			signal_r[k] = signal_daq[k] + coef*acum[k]

		else:  #raw signal just dropped below threshold
		# but raw signal can be negative for a while and still contribute to the
		# reconstructed signal.

			if pulse_on == 1: #reconstructed signal still on
			# switch the pulse off only when recovered signal drops below threshold

				#slide the MAU, still frozen.
				# keep recovering signal
				MAU[k] = MAU[k-1]
				signal_i[k] =MAU[k-1]
				acum[k] = acum[k-1] + signal_daq[k] - offset
				signal_r[k] = signal_daq[k] + coef*acum[k]

				#if the recovered signal drops before trigger line rec pulse is over!
				if signal_r[k] < trigger_line+thr2:
					wait_over = 1  #start tail compensation
					pulse_on = 0   #recovered pulse is over

			else:  #recovered signal has droped below trigger line
			#need to compensate the tail to avoid drifting due to erros in
			#baseline calculatoin

				if wait_over == 1: #compensating pulse
					# recovered signal and raw signal must be equal within a threshold
					# otherwise keep compensating pluse

					if signal_daq[k-1] < signal_r[k-1] - thr_tr:
						# raw signal still below recovered signal

						# is the recovered signal near offset?
						upper = offset + thr2
						lower = offset - thr2

						if signal_r[k-1] > lower and signal_r[k-1] < upper:
							# we are near offset, activate MAU. signal_i follows rec signal

							signal_i[k] = signal_r[k-1]
							MAU[k] = np.sum(signal_i[k-nm:k])/nm

						else:
							# rec signal not near offset MAU frozen

							MAU[k] = MAU[k-1]
							signal_i[k] = MAU[k-1]

						# keep adding recovered signal until
						# it raises above the raw signal

						acum[k] = acum[k-1] + signal_daq[k] - MAU[k]
						signal_r[k] = signal_daq[k] + coef*acum[k]

					else:  # input signal above recovered signal: we are done
						wait_over = 0
						acum[k] = MAU[k-1]
						signal_r[k] = signal_daq[k]
						signal_i[k] = signal_r[k]
						MAU[k] = np.sum(signal_i[k-nm:k])/nm
				else:
					acum[k] = MAU[k-1]
					signal_r[k] = signal_daq[k]
					signal_i[k] = signal_r[k]
					MAU[k] = np.sum(signal_i[k-nm:k])/nm

	# signal recovered, time to cash

	# time_bin = FP.time_DAQ/FP.time_step
	# Q_i = FP.pulse_area_positive(signal_PE_adc)
	# Q_r = time_bin*FP.pulse_area_threshold(signal_r, 3*FP.NOISE_FEE_adc)

	# amax_daq = np.amax(signal_daq)
	# amax_r = np.amax(signal_r)

	# amin_daq = np.amin(signal_daq)
	# amin_r = np.amin(signal_r)

	# range_r = amax_r - amin_r
	# range_daq = amax_daq - amin_daq

	# print """
	# 	Deconv Simple DAQ
	# 	Q (input pulse) = %7.5g
	# 	Q (recovered pulse) = %7.5g (3*noise_rms cut)
	# 	rms = %7.5g

	# 	range DAQ = %7.2f
	# 	range rec = %7.2f
	# 	Ra =  %7.2f
	# 	"""%(Q_i,Q_r, abs(Q_i -Q_r)/Q_i, range_daq, range_r, range_daq/range_r)

	# El /1.5 es pq Juanjo utiliza 1.5 como marcador de senyales

	print BASELINE
	energia = np.sum(((pulse_f-0.5))*(signal_r-BASELINE))/SPE


	if plot_MAU:
#		plot_signal(signal_t_daq/ns,MAU,
#                title = 'MAU',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='')
		plt.plot(MAU)

	if plot_pulse_f:
#		plot_signal(signal_t_daq/ns,pulse_f,
#                title = 'pulse_f',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='')
		plt.plot((pulse_f-0.5)*1000)

	if plot_tail_f:
#		plot_signal(signal_t_daq/ns,tail_f,
#                title = 'tail_f',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='')
		plt.plot(tail_f)

	if plot_trigger_f:
#		plot_signal(signal_t_daq/ns,trigger_f,
#                title = 'trigger_f',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='')
		plt.plot(trigger_f)

	if plot_signal_i:
#		plot_signal(signal_t_daq/ns,signal_i,
#                title = 'signal_i',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='')
		plt.plot(signal_i)

	if plot_acum:
#		plot_signal(signal_t_daq/ns,acum,
#                title = 'acum',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='')
		plt.plot(acum)

	if plot_signal_daq:
#		plot_signal(signal_t_daq/ns,signal_daq,
#                title = 'signal_daq',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='adc counts')
		plt.plot(signal_daq)

	if plot_signal_r:
#		plot_signal(signal_t_daq/ns,signal_r,
#                title = 'signal_r',
#                signal_start=0, signal_end=len(signal_daq)*FP.time_bin*ns,
#                units='adc counts')
		plt.plot(signal_r-BASELINE)

	return  signal_r,energia


def main():
	path = 'D:/DATOS_DAC/2055/pmt_0_trace_evt_1.txt'
	g=pd.read_csv(path)
	f=g.values
	recons, energia_aux = MauDeconv([(1)], 4096-f, 1.636E-3, \
						thr=150, MAU_WindowSize=256, SPE=21.26)
	print energia_aux

if __name__ == "__main__":
	main()



