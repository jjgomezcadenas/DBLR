"""
Controls parameters
"""

cython_dblr = True

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

iset = 6  
frst_evt = 0
lst_evt = 100
n_sigma=1
path = "/Users/jjgomezcadenas/Documents/Development/NEXT/data/1ns/"
histoPath='/Users/jjgomezcadenas/Documents/Development/NEXT/WORK/EP/Signals/'

CPAR = {}
CPAR['saveHistos'] = True
CPAR['signal_analysis'] = True
CPAR['energy_analysis'] = True
CPAR['print_energy_vectors'] = False

CPAR['Plots'] = False
CPAR['Histograms'] = True
CPAR['EnergyHistograms'] = True

CPAR['plot_BLR'] = False
CPAR['plot_I'] = False
CPAR['plot_V'] = True
CPAR['plot_ADC'] = True
CPAR['plot_PE'] = True
CPAR['plot_spe'] = False
CPAR['plot_spe_fee'] = False
CPAR['plot_s1_mc'] = False
CPAR['plot_s2_mc'] = True
CPAR['plot_s1_pmt'] = False
CPAR['plot_s2_pmt'] = False
CPAR['plot_signal_pmt'] = False
CPAR['plot_s1_fee'] = False
CPAR['plot_s2_fee'] = False
CPAR['plot_signal_fee'] = True
CPAR['plot_R'] = True
CPAR['plot_signal_inv_daq'] = False