"""
A utility module for plots with matplotlib
"""
import matplotlib.pyplot as plt
from system_of_units import *
import numpy as np
import os
import sys


def wait():
  """
  Convenience rename of raw_input
  """
  raw_input("Press a key...")


def HSimple1(x,nbins,title='hsimple',xlabel = '', ylabel = 'Frequency', 
             save=False,filename='hsimple.png', filepath='./'):
  plt.hist(x, nbins, histtype='bar', rwidth=0.6)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  
  if save:
    pathfile = filepath+filename
    print "saving histogram %s in %s"%(filename, pathfile)
    plt.savefig(pathfile, bbox_inches='tight')
    plt.close()
  else:
    plt.figure()


def plot_signal(signal_t,signal, 
                title = 'signal', signal_start=0, signal_end=1e+4, units=''):

  ax1 = plt.subplot(1,1,1)
  ax1.set_xlim([signal_start, signal_end])
  SetPlotLabels(xlabel='t (ns)', ylabel='signal (%s)'%units)
  plt.title(title)
  plt.plot(signal_t, signal)
  plt.figure()

def plot_signal2(signal_t,signal, 
                title = 'signal', signal_start=0, signal_end=1e+7,
                signal_min=-1e+7, signal_max=1e+7, 
                units=''):

  ax1 = plt.subplot(1,1,1)
  ax1.set_xlim([signal_start, signal_end])
  ax1.set_ylim([signal_min, signal_max])
  SetPlotLabels(xlabel='t (ns)', ylabel='signal (%s)'%units)
  plt.title(title)
  plt.plot(signal_t, signal)
  plt.figure()

def pulse_plot(pulse_time, pulse_volt):
  """
  Plots pulse
  """
  plt.plot(pulse_time, pulse_volt)
  plt.show()

def SetPlotLabels(xlabel="", ylabel="",grid=True):
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if grid == True:
    plt.grid(which='both', axis='both')


def PrintSignal(title,signal_t,signal):

  print "%s len: signal_t =%d, signal = %d "%(title,len(signal_t), len(signal))
  #print "signal_t =", signal_t
  #print "signal =", signal
  
def PlotPE(cbin,cnt,len_signal_max, s1_l,s1_r,s2_l,s2_r):
    ax1 = plt.subplot(3,1,1)
    ax1.set_xlim([0, len_signal_max/mus])
    SetPlotLabels(xlabel='t (mus)', ylabel='signal (PES)')
    
    plt.plot(cbin/mus, cnt)

    ax2 = plt.subplot(3,1,2)
    ax2.set_xlim([s1_l/mus, s2_r/mus])
    SetPlotLabels(xlabel='t (mus)', ylabel='signal (PES)')
    
    plt.plot(cbin/mus, cnt)
   
    ax3 = plt.subplot(3,1,3)
    ax3.set_xlim([s2_l, s2_r])
    SetPlotLabels(xlabel='t (mus)', ylabel='signal (PES)')
    
    plt.plot(cbin/ns, cnt)
    plt.show()

def PlotSignals(DSGN,CPLOT):
    """
    Plot signals
    """

    #spe
    if CPLOT['plot_spe'] == True:
        if CPLOT['plot_I']:
            print DSGN['spe_I']
            DSGN['spe_I'].PlotSignal(positive = False, title= "SPE after PMT" )
        
        if CPLOT['plot_V']:
            print DSGN['spe_V']
            DSGN['spe_V'].PlotSignal(positive = False, title= "SPE (V) after PMT " )
        
        if CPLOT['plot_ADC']:
            print DSGN['spe_ADC']
            DSGN['spe_ADC'].PlotSignal(positive = False, title= "SPE (ADC) after PMT" )
    
    if CPLOT['plot_spe_fee']:


        print DSGN['spe_FEE']
        DSGN['spe_FEE_NN'].PlotSignal(positive = False, title= "SPE after FEE NN" )
        DSGN['spe_FEE'].PlotSignal(positive = False, title= "SPE after FEE" )

        print DSGN['spe_DAQ']
        DSGN['spe_DAQ'].PlotSignal(positive = False, title= "SPE after DAQ" )


    # s1 and s2 MC average of 12 channels
    if CPLOT['plot_s1_mc'] == True:
        print DSGN['s1_MC']
        DSGN['s1_MC'].PlotSignal(positive = False, title= "S1 MC") 
        DSGN['s1_MC'].PlotSignal(positive = True, title= "S1 MC: positive downscaled") 
        

    if CPLOT['plot_s2_mc'] == True:
        print DSGN['s2_MC']
        DSGN['s2_MC'].PlotSignal(positive = False, title= "S2 MC") 
        DSGN['s2_MC'].PlotSignal(positive = True, title= "S2 MC: positive downscaled") 
       

    #s1 and s2 signal from PMT
    if CPLOT['plot_s1_pmt'] == True:
        if CPLOT['plot_I']:
            print DSGN['s1_PMT']
            DSGN['s1_PMT'].PlotSignal(positive = False, title= "S1 after PMT" )
            DSGN['s1_PMT'].PlotSignal(positive = True, title= "S1 after PMT: pos downscaled" )
        
        if CPLOT['plot_V']:
            print DSGN['s1_PMT_V']
            DSGN['s1_PMT_V'].PlotSignal(positive = False, title= "S1 (V) after PMT" )
            DSGN['s1_PMT_V'].PlotSignal(positive = True, title= "S1 (V) after PMT: pos downscaled" )
        
        if CPLOT['plot_ADC']:
            print DSGN['s1_PMT_ADC']
            DSGN['s1_PMT_ADC'].PlotSignal(positive = False, title= "S1 (ADC) after PMT" )
            DSGN['s1_PMT_ADC'].PlotSignal(positive = True, title= "S1 (ADC) after PMT: pos downscaled")

    if CPLOT['plot_s2_pmt'] == True:
        if CPLOT['plot_I']:
            print DSGN['s2_PMT']
            DSGN['s2_PMT'].PlotSignal(positive = False, title= "S2 after PMT" )
            DSGN['s2_PMT'].PlotSignal(positive = True, title= "S2 after PMT: pos downscaled")
       
        if CPLOT['plot_V']:
            print DSGN['s2_PMT_V']
            DSGN['s2_PMT_V'].PlotSignal(positive = False, title= "S2 (V) after PMT" )
            DSGN['s2_PMT_V'].PlotSignal(positive = True, title= "S2 (V) after PMT: pos downscaled")
        
        if CPLOT['plot_ADC']:
            print DSGN['s2_PMT_ADC']
            DSGN['s2_PMT_ADC'].PlotSignal(positive = False, title= "S2 (ADC) after PMT" )
            DSGN['s2_PMT_ADC'].PlotSignal(positive = True, title= "S2 (ADC) after PMT: pos downscaled")
        
    if CPLOT['plot_signal_pmt'] == True:

        print DSGN['signal_PMT_ADC']
        DSGN['signal_PMT_ADC'].PlotSignal(positive = False, title= "Signal (ADC) after PMT" )
        DSGN['signal_PMT_ADC'].PlotSignal(positive = True, 
                                          title= "Signal (ADC) after PMT: pos downscaled")
        
    #s1 and s2 signal from FEE
    if CPLOT['plot_s1_fee'] == True:
        print DSGN['s1_FEE_NN']
        DSGN['s1_FEE_NN'].PlotSignal(positive = False, title= "S1 after FEE NN" )
        
        print DSGN['s1_FEE']
        DSGN['s1_FEE'].PlotSignal(positive = False, title= "S1 after FEE " )
        
        print DSGN['s1_DAQ']
        DSGN['s1_DAQ'].PlotSignal(positive = False, title= "S1 after DAQ " )

    if CPLOT['plot_s2_fee'] == True:
        print DSGN['s2_FEE']
        DSGN['s2_FEE'].PlotSignal(positive = False, title= "S2 after FEE" )
        
        print DSGN['s2_DAQ']
        DSGN['s2_DAQ'].PlotSignal(positive = False, title= "S2 after DAQ" )

    if CPLOT['plot_signal_fee'] == True:
        print DSGN['signal_DAQ']
        DSGN['signal_DAQ'].PlotSignal(positive = False, title= "Signal after DAQ" )
        print DSGN['signal_DAQ_off']
        DSGN['signal_DAQ_off'].PlotSignal(positive = False, title= "Signal after DAQ offset" )
        
    #recovered
    if CPLOT['plot_R'] == True:
        print DSGN['signal_R']
        DSGN['signal_R'].PlotSignal(positive = False, title= "Signal Recovered" )
        
        # DSGN['s2_R'].PlotSignal(positive = False, title= "S2 Recovered" )
        # print DSGN['s2_R']

