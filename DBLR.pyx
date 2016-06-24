import numpy as np
from scipy import signal as SGN
import FEParam as FP

import sys
from system_of_units import *

cdef int MAU_WindowSize = FP.MAU_WindowSize
cdef float time_DAQ = FP.time_bin



def BLR(double[:] signal_t_daq, double[:] signal_daq, double coef, 
        int n_sigma = 3, double NOISE_ADC=0.5):
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

    cdef int len_signal_daq = len(signal_daq)
    MAU = np.zeros(len_signal_daq, dtype=np.double)
    acum = np.zeros(len_signal_daq, dtype=np.double)
    signal_r = np.zeros(len_signal_daq, dtype=np.double)
    
    signal_i = np.copy(signal_daq) #uses to update MAU while procesing signal

    cdef double thr = n_sigma*NOISE_ADC
    cdef double thr_cmp = 2*thr  # to follow baseline tail
    cdef double thr_tr = thr/5. # to conclude BLR when signal_deconv = signal_raw
    
    #MAU_WindowSize = 40 # provisional
    cdef int nm = MAU_WindowSize
    B_MAU       =   (1./nm)*np.ones(nm)

#   MAU averages the signal in the initial tranch 
#    allows to compute the baseline of the signal  
    
    MAU[0:nm] = SGN.lfilter(B_MAU,1, signal_daq[0:nm])
    acum[nm] =  MAU[nm]

#----------

# While MAU inits BLR is switched off, thus signal_r = signal_daq 

    signal_r[0:nm] = signal_daq[0:nm] 
    cdef int pulse_on=0
    cdef int wait_over=0
    cdef double offset = 0

    # MAU has computed the offset using nm samples
    # now loop until the end of DAQ window
    cdef int k
    for k in range(nm,len_signal_daq): 

        trigger_line = MAU[k-1] + thr

        # condition: raw signal raises above trigger line and 
        # we are not in the tail
        # (wait_over == 0)
        if signal_daq[k] > trigger_line and wait_over == 0:

            # if the pulse just started pulse_on = 0.
            # In this case compute the offset as value 
            #of the MAU before pulse starts (at k-1)

            if pulse_on == 0: # pulse just started
                #offset computed as the value of MAU before pulse starts
                offset = MAU[k-1]  
                pulse_on = 1 

            #Freeze the MAU
            MAU[k] = MAU[k-1]  
            signal_i[k] =MAU[k-1]  #signal_i follows the MAU
            
            #update recovered signal, correcting by offset
            acum[k] = acum[k-1] + signal_daq[k] - offset;
            signal_r[k] = signal_daq[k] + coef*acum[k] 
            
        else:  #raw signal just dropped below threshold
        # but raw signal can be negative for a while and still contribute to the
        # reconstructed signal.

            if pulse_on == 1: #reconstructed signal still on
            # switch the pulse off only when recovered signal 
            #drops below threshold
                
                #slide the MAU, still frozen. 
                # keep recovering signal
                MAU[k] = MAU[k-1]  
                signal_i[k] =MAU[k-1]
                acum[k] = acum[k-1] + signal_daq[k] - offset;
                signal_r[k] = signal_daq[k] + coef*acum[k] 
            
                #if the recovered signal drops before trigger line 
                #rec pulse is over!
                if signal_r[k] < trigger_line:
                    wait_over = 1  #start tail compensation
                    pulse_on = 0   #recovered pulse is over

            else:  #recovered signal has droped below trigger line
            #need to compensate the tail to avoid drifting due to erros in 
            #baseline calculatoin

                if wait_over == 1: #compensating pulse
                    # recovered signal and raw signal 
                    #must be equal within a threshold
                    # otherwise keep compensating pluse


                    if signal_daq[k-1] < signal_r[k-1] - thr_tr:
                        # raw signal still below recovered signal 

                        # is the recovered signal near offset?
                        upper = offset + thr_cmp
                        lower = offset - thr_cmp

                        if signal_r[k-1] > lower and signal_r[k-1] < upper:
                            # we are near offset, activate MAU. 
                            #signal_i follows rec signal
                            
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
                    wait_over = 0
                    acum[k] = MAU[k-1]
                    signal_r[k] = signal_daq[k]
                    signal_i[k] = signal_r[k]
                    MAU[k] = np.sum(signal_i[k-nm:k])/nm
                        
    return  signal_r

def FindSignalAboveThr(double[:] signal_t, double[:] signal, double threshold = 0.):
    """
    Finds positive signals above threshold. 
    """    
   
    cdef double pulse_on = 0
    cdef int len_signal = len(signal)
    cdef int k

    pulse_f = np.zeros(len_signal, dtype=np.double)
    t_f = np.zeros(len_signal, dtype=np.double)

    cdef int i=0
    for k in range(len_signal): 

        if signal[k] > threshold and pulse_on == 0:
            pulse_f[i] = signal[k]
            t_f[i] = signal_t[k]
            pulse_on = 1
            i+=1
        elif signal[k] < threshold and pulse_on == 1:
            pulse_f[i] = 0.
            t_f[i] =signal_t[k]
            pulse_on = 0
            i+=1
        elif signal[k] > threshold and pulse_on == 1:
            pulse_f[i] = signal[k]
            t_f[i] = signal_t[k]
            i+=1
        
    
    signal_f = np.zeros(i, dtype=np.double)
    time_f = np.zeros(i, dtype=np.double)

    signal_f[:] = pulse_f[0:i]
    time_f[:] = t_f[0:i]
    
    return time_f, signal_f
