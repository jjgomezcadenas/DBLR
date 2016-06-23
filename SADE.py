"""
SADE: Signal Analysis and DEconvolution
Pure python version
JJGC June 2016
"""

from math import *
from PlotUtil import *
from SadeSignals import *

import numpy as np
import os
import sys

import FEParam as FP
import SPE as SP
import FEE as FE

from scipy import signal as SGN
import tables as tb

len_signal_max = 1200000


dataset_name = ["Kr83.z5cm.1ns","electrons.511keV.1ns","electrons.1275keV.1ns",
                "electrons.2615keV.1ns","bb0nu.1ns","electrons.1275keV.1ns.15bar",
                "electrons.2615keV.1ns.15bar","electrons.511keV.1ns.15bar"]

dataset_s1l =[599000,599000,595000,595000,595000,595000,595000,595000]
dataset_s1r =[600500,600500,600500,600500,600500,600500,600500,600500]
dataset_s2l =[600500,601000,601000,601000,601000,601000,601000,601000]
dataset_s2r =[680000,720000,1000000,1100000,1120000,1100000,1100000,720000]

saveHistos = True
CPLOT={}
CPLOT['Plots'] = False
CPLOT['Areas'] = False
CPLOT['Histograms'] = True

CPLOT['plot_I'] = False
CPLOT['plot_V'] = True
CPLOT['plot_ADC'] = True
CPLOT['plot_PE'] = True
CPLOT['plot_spe'] = False
CPLOT['plot_spe_fee'] = False
CPLOT['plot_s1_mc'] = False
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
    
    ds = SignalDef(dataset_name[iset], dataset_s1l[iset], dataset_s1r[iset], 
                   dataset_s2l[iset], dataset_s2r[iset])
    #ds = SignalDef("Kr83.z5cm.1ns", 599000, 600500, 600500, 680000)
    #ds = SignalDef("electrons.511keV.1ns", 599000, 600500, 601000, 720000)
    #ds = SignalDef("electrons.1275keV.1ns", 595000, 600500, 610000, 1000000)
    #ds = SignalDef("electrons.2615keV.1ns", 595000, 600500, 610000, 1100000)

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
    histo_t = np.arange(0.0, len_signal_max*1.+0.0, 1., dtype='Float32')

    #reduced histogram for signal (clip left part of signal, include s1 and s2)
    signal_t = np.arange(0.0, ds.len_signal*1.+0.0, 1., dtype='Float32')

    #reduced histogram for s2 (clip left part of signal, from zero to t_S1 included)
    s2_t = np.arange(0.0, ds.len_s2*1.+0.0, 1., dtype='Float32')

    #reduced histogram for S1 (clip left and right part of signal, leaving a window for S1)
    s1_t = np.arange(0.0, ds.len_s1*1.+0.0, 1., dtype='Float32')
    
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
   
        signal_dec = MauDeconv(signal_t_daq, signal_daq, coef, n_sigma = 2)
        #s2_dec = MD.MauDeconv(s2t_daq, s2_daq, coef, n_sigma = 2)

        DSGN['spe_I'] = Signal(t_spe, pulse_spe,signal_name ='SPE', threshold = 0.1*muA, 
                        units='muA', downscale=False)
        DSGN['spe_V'] = Signal(t_spe, signal_spe,signal_name ='SPE',threshold = 0.1*mV, 
                        units='mV', downscale=False)
        DSGN['spe_ADC'] = Signal(t_spe, signal_spe/FP.voltsToAdc,signal_name ='SPE',
                                threshold = 0.1*adc, units='adc', downscale=False)

        DSGN['spe_FEE'] = Signal(t_spe, spe_fee, signal_name ='SPE, FEE response (V)',
                                    threshold = 0, units='mV', downscale=False)

        DSGN['spe_FEE_NN'] = Signal(t_spe, spe_fee_no_noise, 
                                    signal_name ='SPE, FEE response (no noise) (V)',
                                    threshold = 0, units='mV', downscale=False)

        DSGN['spe_DAQ'] = Signal(spe_t_daq, spe_DAQ, signal_name ='SPE, FEE response (ADC)',
                                threshold = 0., units='adc', stype='DAC', downscale=False)

        DSGN['s2_MC'] = Signal(s2_t, s2_histo_PE, signal_name ='S2 12 chan average (PES)',
                               threshold = 0.01*pes, units = 'pes')
        DSGN['s1_MC'] = Signal(s1_t, s1_histo_PE, signal_name ='S1 12 chan average (PES)',
                               threshold = 0.001*pes, units = 'pes', downscale=False)
        DSGN['s2_PMT'] = Signal(s2_t, s2_PE,signal_name ='S2, PMT response (I)',
                                threshold = 0.01*muA, units='muA')
        DSGN['s1_PMT'] = Signal(s1_t, s1_PE, signal_name ='S1, PMT response (I)', 
                                units='muA', downscale=False)
        DSGN['s2_PMT_V'] = Signal(s2_t, fee.VSignal(s2_PE),signal_name ='S2, PMT response (V)',
                                  threshold = 0.01*mV, units='mV')
        DSGN['s1_PMT_V'] = Signal(s1_t, fee.VSignal(s1_PE), signal_name ='S1, PMT response (V)',
                                  units='mV', downscale=False)

        DSGN['signal_PMT_ADC'] = Signal(signal_t, fee.VSignal(signal_PE)/FP.voltsToAdc,
                                    signal_name ='Signal, PMT response (ADC)',
                                    threshold = 0.01*adc, units='adc')

        DSGN['s2_PMT_ADC'] = Signal(s2_t, fee.VSignal(s2_PE)/FP.voltsToAdc,
                                    signal_name ='S2, PMT response (ADC)',
                                    threshold = 0.01*adc, units='adc')
        DSGN['s1_PMT_ADC'] = Signal(s1_t, fee.VSignal(s1_PE)/FP.voltsToAdc, 
                                    signal_name ='S1, PMT response (ADC)',
                                    threshold = 0.01*adc,
                                    units='adc', downscale=False)

        DSGN['s2_FEE'] = Signal(s2_t, s2_fee,threshold = 0, 
                                signal_name ='S2, FEE response (V)',
                                units='mV')
        DSGN['s1_FEE_NN'] = Signal(s1_t, s1_fee_no_noise,
                                    signal_name ='S1, FEE response, no noise (V)',
                                    threshold = 0, units='mV', downscale=False)
        DSGN['s1_FEE'] = Signal(s1_t, s1_fee, signal_name ='S1, FEE response (V)',
                                    threshold = 0, units='mV', downscale=False)

        DSGN['signal_DAQ_off'] = Signal(signal_t_daq, signal_daq -FP.offset, 
                                    signal_name ='Signal, DAQ offset',
                                    threshold = 0.0, units='adc', stype='DAC', downscale=False)

        DSGN['s2_DAQ'] = Signal(s2t_daq, s2_daq, signal_name ='S2, FEE response (ADC)',
                                threshold = 0., units='adc', stype='DAC', downscale=False)
        DSGN['s1_DAQ'] = Signal(s1t_daq, s1_daq, signal_name ='S1, FEE response (ADC)',
                                threshold = 0., units='adc', stype='DAC', downscale=False)

        DSGN['signal_DAQ'] = Signal(signal_t_daq, signal_daq, 
                                signal_name ='Signal, DAQ (ADC)',
                                threshold = 0.1, units='adc', stype='DAC', downscale=False)
        DSGN['signal_R'] = Signal(signal_t_daq, signal_dec, 
                                signal_name ='Signal, Deconvoluted (ADC)',
                                threshold = 0.1, units='adc', stype='DAC', downscale=False)

        # DSGN['s2_R'] = Signal(s2t_daq, s2_dec, signal_name ='S2, Deconvoluted (ADC)',
        #                         threshold = 0., units='adc', stype='DAC', downscale=False)

        # DSGN['s1_R'] = Signal(s1t_daq, s1_dec, signal_name ='S1, Deconvoluted (ADC)',
        #                         threshold = 0., units='adc', downscale=False)


        if CPLOT['Plots']:
            PlotSignals(DSGN,CPLOT)

        if CPLOT['Areas']:
            ComputeAreas(DSGN)
        
        CSGN.append(DSGN)
        

    SignalAnalysis(CSGN,CPLOT,saveHistos, filepath =hpath)

    if saveHistos == False or CPLOT['Histograms'] == False:
        plt.show()
    #plt.draw()
    


if __name__ == '__main__':

    import cProfile
    iset = 6  #iset: Kr =0, e511=1, e1275=2, e2651=3, bb0nu=4
    frst_evt = 0
    lst_evt = 10
    path = "/Users/jjgomezcadenas/Documents/Development/NEXT/data/1ns/"
    histoPath='/Users/jjgomezcadenas/Documents/Development/NEXT/WORK/EP/Signals/'
    SADE(path,histoPath,iset,nmin=frst_evt,nmax=lst_evt)
    #cProfile.run('SADE(iset,nmax=maxevt)', sort='time')

        
     