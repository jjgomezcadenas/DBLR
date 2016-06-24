"""
SADE: Signal Analysis and DEconvolution
JJGC June 2016
"""

from PlotUtil import *
from SadeSignals import *
import DBLR as DB

import numpy as np
import os
import sys

import FEParam as FP
import SPE as SP
import FEE as FE

from scipy import signal as SGN
import tables as tb

len_signal_max = 1200000


dataset_name = ["Kr83.z5cm.1ns","Kr83.z10cm.1ns","Kr83.z15cm.1ns",
                "electrons.511keV.1ns","electrons.1275keV.1ns",
                "electrons.2615keV.1ns",
                "bb0nu.1ns",
                "Kr83.1ns.15bar",
                "electrons.1275keV.1ns.15bar",
                "electrons.2615keV.1ns.15bar","electrons.511keV.1ns.15bar",
                "bb0nu.1ns.15bar"]

dataset_s1l = 599000
dataset_s1r =600500 
dataset_s2l =600500
dataset_s2r =1100000

saveHistos = True
CPLOT={}
CPLOT['Plots'] = True
CPLOT['Areas'] = True
CPLOT['Histograms'] = False

CPLOT['plot_I'] = False
CPLOT['plot_V'] = True
CPLOT['plot_ADC'] = True
CPLOT['plot_PE'] = True
CPLOT['plot_spe'] = True
CPLOT['plot_spe_fee'] = True
CPLOT['plot_s1_mc'] = True
CPLOT['plot_s2_mc'] = True
CPLOT['plot_s1_pmt'] = False
CPLOT['plot_s2_pmt'] = False
CPLOT['plot_signal_pmt'] = False
CPLOT['plot_s1_fee'] = False
CPLOT['plot_s2_fee'] = False
CPLOT['plot_signal_fee'] = True
CPLOT['plot_R'] = True
CPLOT['plot_signal_inv_daq'] = False

def SADE(path,histoPath,iset,nmin=0,nmax=100):
    
    ds = SignalDef(dataset_name[iset], dataset_s1l, dataset_s1r, 
                   dataset_s2l, dataset_s2r)
    
    hdir = histoPath+dataset_name[iset]
    if not os.path.exists(hdir):
        os.makedirs(hdir)
    hpath = hdir+'/'

    print """
    Reading h5 file  = %s
    """%(ds.h5file)

    print "histograms saved in directory %s"%hdir

    FP.print_FEE()
    wait()


    #time vector

    #full histogram (1.2E+6 bins)
    histo_t = np.arange(0.0, len_signal_max*1.+0.0, 1., dtype=np.double)

    #reduced histogram for signal (clip left part of signal, include s1 and s2)
    signal_t = np.arange(0.0, ds.len_signal*1.+0.0, 1., dtype=np.double)

    #reduced histogram for s2 (clip left part of signal, from zero to t_S1 included)
    s2_t = np.arange(0.0, ds.len_s2*1.+0.0, 1., dtype=np.double)

    #reduced histogram for S1 (clip left and right part of signal, leaving a window for S1)
    s1_t = np.arange(0.0, ds.len_s1*1.+0.0, 1., dtype=np.double)
    
    s2_histo_PE = np.zeros(ds.len_s2)  #number of PE for S2
    s1_histo_PE = np.zeros(ds.len_s1)  # number of PE for S1
    signal_histo_PE = np.zeros(ds.len_signal)  #number of PE for S2

    #objects defining the electronics
    spe = SP.SPE(pmt_gain=FP.PMT_GAIN,x_slope = 5*ns,x_flat = 1*ns)
    fee = FE.FEE(PMTG=FP.PMT_GAIN, C=FP.C,R= FP.R, 
                 f=FP.freq_LPF, fn=FP.freq_HPF, RG=FP.V_GAIN)
    daq = FE.DAQ(NBITS=FP.NBITS, NBITS_FRAC=FP.NBITS_FRAC,
                 time_sample=FP.time_DAQ, 
                 LSB = 2*volt)

    #open h5file
    h5file = tb.open_file(path+ds.h5file, mode = "r")
    print(h5file)
    
    i=0
    
    CSGN = [] #container of signals 

    for array in h5file.walk_nodes("/", "CArray"):
        i+=1
        if i > nmax:
            break 

        DSGN ={} #dictionary of signals
        
        print " analyzing event  = %d"%i
        # print(array)

        cnt = array.read()  #load histogram from h5 file

        if i <= nmin:
            continue
        
        if CPLOT['Plots'] and CPLOT['plot_PE']:
            PlotPE(histo_t,cnt,len_signal_max,ds.s1_l,ds.s1_r,ds.s2_l,ds.s2_r)

        #s2 and s1 histograms from True MC (average of PMTs, thus fractions allowed) in PES 
        s2_histo_PE[:] = cnt[ds.s2_l:ds.s2_r]
        s1_histo_PE[:] = cnt[ds.s1_l:ds.s1_r]
        signal_histo_PE[:] = cnt[ds.sgn_l:ds.sgn_r]

        #s2 and s1 signals
        s2_PE = spe.SpePulseFromVectorPE(s2_histo_PE) #PMT response
        s1_PE = spe.SpePulseFromVectorPE(s1_histo_PE) #PMT response
        signal_PE = spe.SpePulseFromVectorPE(signal_histo_PE) #PMT response

        #single photoelectron (spe)
        t0 = len(s1_t)*ns/2
        tmax = len(s1_t)*ns
        t_spe, pulse_spe = spe.SpePulse(t0, tmax = tmax)
        signal_spe = fee.VSignal(pulse_spe)
        
    
        spe_fee_no_noise = fee.FEESignal(pulse_spe, noise_rms=0)
        # print "area of spe = %7.2f"%(np.sum(spe_fee_no_noise)/mV)
        # wait()

        spe_fee = fee.FEESignal(pulse_spe, noise_rms=FP.NOISE_FEE_rms) #in single PMT
        spe_t_daq,spe_DAQ = daq.DAQSignal(t_spe, spe_fee, noise_rms=0)

        # effect of electronics: pass filters and output voltage
        # add noise
        signal_fee = fee.FEESignal(signal_PE, noise_rms=FP.NOISE_FEE) #in average
        s2_fee = fee.FEESignal(s2_PE, noise_rms=FP.NOISE_FEE) #in average
        s1_fee_no_noise = fee.FEESignal(s1_PE, noise_rms=0)
        s1_fee = fee.FEESignal(s1_PE, noise_rms=FP.NOISE_FEE_rms) #in single PMT
        
        #Signal out of DAQ
        signal_t_daq, signal_daq = daq.DAQSignal(signal_t, signal_fee, noise_rms=0) 
        s2t_daq, s2_daq = daq.DAQSignal(s2_t, s2_fee, noise_rms=0) #nose added to s2_fee
        s1t_daq, s1_daq = daq.DAQSignal(s1_t, s1_fee, noise_rms=0) #nose added to s1_fee
        
        #Deconvolution
        signal_inv_daq = fee.InverseSignalDAQ(signal_t)  #inverse function
        coef = signal_inv_daq[10]  #accumulator coefficient
        #print "inverse coef fee: = %7.2g "%coef

        if CPLOT['plot_signal_inv_daq'] and CPLOT['Plots'] :   
            plot_signal(signal_t/ns,signal_inv_daq,
                title = 'Inverse DAQ', 
                signal_start=0*ns, signal_end=10*ns, 
                units='')
            print "inverse coef fee: = %7.2g "%coef

        # print "calling MauDeconv"
   
        #signal_dec = MauDeconv(signal_t_daq, signal_daq, coef, n_sigma = 2)
        #s2_dec = MD.MauDeconv(s2t_daq, s2_daq, coef, n_sigma = 2)
        signal_dec = DB.BLR(signal_t_daq, signal_daq, coef, n_sigma = 2, NOISE_ADC=0.5)

        DSGN['spe_I'] = Signal(t_spe, pulse_spe, threshold = 0.01*muA)
        DSGN['spe_V'] = Signal(t_spe, signal_spe,threshold = 0.01*mV)
        DSGN['spe_ADC'] = Signal(t_spe, signal_spe/FP.voltsToAdc, threshold = 0.01*adc)

        DSGN['spe_FEE'] = Signal(t_spe, spe_fee, threshold = 0)
        DSGN['spe_FEE_NN'] = Signal(t_spe, spe_fee_no_noise, threshold = 0)
        DSGN['spe_DAQ'] = Signal(spe_t_daq, spe_DAQ, threshold = 0,  stype='DAC')

        DSGN['s2_MC'] = Signal(s2_t, s2_histo_PE,threshold = 0.01*pes, downscale=True)
        DSGN['s1_MC'] = Signal(s1_t, s1_histo_PE, threshold = 0.001*pes, downscale=True)
                               
        DSGN['s2_PMT'] = Signal(s2_t, s2_PE,threshold = 0.01*muA, downscale=True)
        DSGN['s1_PMT'] = Signal(s1_t, s1_PE, threshold = 0.01*muA, downscale=True)
        DSGN['s2_PMT_V'] = Signal(s2_t, fee.VSignal(s2_PE),threshold = 0.01*mV, downscale=True)
        DSGN['s1_PMT_V'] = Signal(s1_t, fee.VSignal(s1_PE), threshold = 0.01*mV, downscale=True)

        DSGN['signal_PMT_ADC'] = Signal(signal_t, fee.VSignal(signal_PE)/FP.voltsToAdc,
                                    threshold = 0.01*adc, downscale=True)

        DSGN['s2_PMT_ADC'] = Signal(s2_t, fee.VSignal(s2_PE)/FP.voltsToAdc,
                                    threshold = 0.01*adc, downscale=True)
        DSGN['s1_PMT_ADC'] = Signal(s1_t, fee.VSignal(s1_PE)/FP.voltsToAdc, 
                                    threshold = 0.01*adc,downscale=False)

        DSGN['s2_FEE'] = Signal(s2_t, s2_fee,threshold = 0, downscale=False)
        DSGN['s1_FEE_NN'] = Signal(s1_t, s1_fee_no_noise,threshold = 0, downscale=False)
        DSGN['s1_FEE'] = Signal(s1_t, s1_fee, threshold = 0, downscale=False)
                                    

        DSGN['signal_DAQ_off'] = Signal(signal_t_daq, signal_daq -FP.offset, 
                                    threshold = 0, downscale=False)

        DSGN['s2_DAQ'] = Signal(s2t_daq, s2_daq, threshold = 0, downscale=False)
        DSGN['s1_DAQ'] = Signal(s1t_daq, s1_daq, 
                                threshold = 0, stype='DAC', downscale=False)

        DSGN['signal_DAQ'] = Signal(signal_t_daq, signal_daq, 
                                threshold = 0, stype='DAC', downscale=False)
        DSGN['signal_R'] = Signal(signal_t_daq, signal_dec, 
                                threshold = 0, stype='DAC', downscale=False)


        if CPLOT['Plots']:
            PlotSignals(DSGN,CPLOT)

        if CPLOT['Areas']:
            ComputeAreas(DSGN)
        
        CSGN.append(DSGN)
        

    SignalAnalysis(CSGN,CPLOT,saveHistos, filepath =hpath)

    if saveHistos == False or CPLOT['Histograms'] == False:
        plt.show()
    


if __name__ == '__main__':

    import cProfile
    iset = 0  
    frst_evt = 0
    lst_evt = 10
    path = "/Users/jjgomezcadenas/Documents/Development/NEXT/data/1ns/"
    histoPath='/Users/jjgomezcadenas/Documents/Development/NEXT/WORK/EP/Signals/'
    SADE(path,histoPath,iset,nmin=frst_evt,nmax=lst_evt)
    #cProfile.run('SADE(iset,nmax=maxevt)', sort='time')

        
     