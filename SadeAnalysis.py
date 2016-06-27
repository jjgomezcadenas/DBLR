"""
Analysis 
"""
from PlotUtil import *
import numpy as np

def EnergyAnalysis(SE,CPAR, filepath ="./"):
    """
    Takes an instance of Sade Energy, computes histograms
    """

    print """
        Sade Energy = %s
    """%(SE)

    epes25ns = SE.RecoveredEnergyInPES(ebin='25ns')
    
    qr25ns = 1000*SE.RatioTrueToRecovered(ebin='25ns')
    qr1ns = 1000*SE.RatioTrueToRecovered(ebin='1ns')

    etrue25ns = SE.TrueEnergyInPES(ebin='25ns')
    etrue1ns = SE.TrueEnergyInPES(ebin='1ns')

    # print """

    #     True Energy (MC only FANO): 25 ns
    #     mean = %7.2f std = %7.2f, sigma_E/E (FWHM) = %7.2f

    #     True Energy (MC only FANO): 1 ns
    #     mean = %7.2f std = %7.2f, sigma_E/E (FWHM) = %7.2f 

    #     Difference Ratio (x1000) = %7.2g

    #     Reconstructed energy (in PES, 25 ns includes electronics noise, recovery)
    #     mean = %7.2f std = %7.2f, sigma_E/E (FWHM) = %7.2f

    #     Difference Ratio True to Recovered (25 ns) x 1000 = %7.2f
    #     Difference Ratio True to Recovered (1 ns) x 1000 = %7.2f

    # """%(len(epes25ns), 
    #      np.average(etrue25ns), np.std(etrue25ns),
    #      ResFWHM(np.std(etrue25ns),np.average(etrue25ns)),
    #      np.average(etrue1ns), np.std(etrue1ns),
    #      ResFWHM(np.std(etrue1ns),np.average(etrue1ns)),
    #      differenceRatio(np.average(etrue1ns),np.average(etrue25ns)),
    #      np.average(epes25ns), np.std(epes25ns),
    #      ResFWHM(np.std(epes25ns),np.average(epes25ns)),
    #      np.average(qr25ns),
    #      np.average(qr1ns)
    #      )

    if CPAR['EnergyHistograms'] == True:

        bins = hbins(etrue25ns, nsigma=5, nbins=20)
        HSimple1(etrue25ns,bins,title="True energy in PES (25 ns)",xlabel = "pes",
            save=CPAR['saveHistos'],filename='etrue25ns.png', filepath=filepath)

        bins = hbins(epes25ns, nsigma=5, nbins=20)
        HSimple1(epes25ns,bins,title="Recovered energy in PES (25 ns)",xlabel = "pes",
            save=CPAR['saveHistos'],filename='epes25ns.png', filepath=filepath)

        bins = hbins(qr25ns, nsigma=5, nbins=20)
        HSimple1(qr25ns,bins,title="Rec ratio (25 ns)",xlabel = "",
            save=CPAR['saveHistos'],filename='qr25ns.png', filepath=filepath)

def SignalAnalysis(CSGN,CPAR, filepath ="./"):
    """
    Histogram main properties
    """
    s1Area = np.zeros(len(CSGN))
    s2Area = np.zeros(len(CSGN))
    s2Length = np.zeros(len(CSGN))
    s2Peak = np.zeros(len(CSGN))
    s2Avg = np.zeros(len(CSGN))

    SPEAdcPeak = np.zeros(len(CSGN))
    SPEAdcArea = np.zeros(len(CSGN))

    SignalRecPeak = np.zeros(len(CSGN))
    SignalRecArea = np.zeros(len(CSGN))
       
    SignalDAQPeak = np.zeros(len(CSGN))
    SignalDAQArea = np.zeros(len(CSGN))

    SignalDAQOffMax = np.zeros(len(CSGN))
    SignalDAQOffMin = np.zeros(len(CSGN))
    

    i = 0
    for DSGN in CSGN:

        SPEAdcPeak[i] = DSGN['spe_ADC'].Peak()
        SPEAdcArea[i] = DSGN['spe_ADC'].Area()

        s1Area[i] = DSGN['s1_MC'].Area()
        s2Area[i] = DSGN['s2_MC'].Area()
        s2Length[i] = DSGN['s2_MC'].Length(units='mus')
        s2Peak[i] = DSGN['s2_MC'].Peak()
        s2Avg[i] = DSGN['s2_MC'].Avg()

        SignalDAQPeak[i] = DSGN['signal_DAQ'].Peak()
        SignalDAQArea[i] = DSGN['signal_DAQ'].Area()

        SignalRecPeak[i] = DSGN['signal_R2'].Peak()
        SignalRecArea[i] = DSGN['signal_R2'].Area()
        
        SignalDAQOffMax[i] = DSGN['signal_DAQ_off'].Peak()
        SignalDAQOffMin[i] = DSGN['signal_DAQ_off'].Minimum()

    
        i+=1 
    
    print """

    SPE DAQ
    peak 25 ns = %7.2f adc counts
    area = %7.2g adc counts

    S1
    area  = %7.2f pes

    S2 
    peak  = %7.2f pes
    average heigth  = %7.2f pes
    area  = %7.2f pes
    Length  = %7.2f mus

    DAQ 
    peak  = %7.2f adc
    area  = %7.2f adc

    DAQ offset
    maximum  = %7.2f adc
    minimum  = %7.2f adc

    REC
    peak  = %7.2f adc
    area  = %7.2f adc

    
    """%(np.average(SPEAdcPeak),np.average(SPEAdcArea),
        np.average(s1Area),
        np.average(s2Peak),np.average(s2Avg),np.average(s2Area),
        np.average(s2Length),
        np.average(SignalDAQPeak), np.average(SignalDAQArea),
        np.average(SignalDAQOffMax),np.average(SignalDAQOffMin),
        np.average(SignalRecPeak),np.average(SignalRecArea)
        )


    if CPAR['Histograms']:

        
        bins = hbins(s1Area, nsigma=5, nbins=20)
        HSimple1(s1Area,bins,title="S1 area in PES (1ns)",xlabel = "pes",
            save=CPAR['saveHistos'],filename='s1Area.png', filepath=filepath)

        bins = hbins(s2Area, nsigma=5, nbins=20)
        HSimple1(s2Area,bins,title="S2 area (PES)",xlabel = "pes",
            save=CPAR['saveHistos'],filename='s2Area.png', filepath=filepath)

        bins = hbins(s2Length, nsigma=5, nbins=20)
        HSimple1(s2Length,bins,title="S2 length (mus)",xlabel = "mus",
            save=CPAR['saveHistos'],filename='s2Length.png', filepath=filepath)

        bins = hbins(s2Peak, nsigma=5, nbins=20)
        HSimple1(s2Peak,bins,title="S2 peak in PES (25 ns)",xlabel = "pes",
            save=CPAR['saveHistos'],filename='s2Peak.png', filepath=filepath)

        bins = hbins(s2Avg, nsigma=5, nbins=20)
        HSimple1(s2Avg,bins,title="S2 avg in PES (25 ns)",xlabel = "pes",
            save=CPAR['saveHistos'],filename='s2Avg.png', filepath=filepath)

        bins = hbins(SignalRecPeak, nsigma=5, nbins=20)
        HSimple1(SignalRecPeak,bins,title="Signal R  peak in ADC (25 ns)",xlabel = "adc",
            save=CPAR['saveHistos'],filename='SignalRecPeak.png', filepath=filepath)
        
        bins = hbins(SignalDAQPeak, nsigma=5, nbins=20)
        HSimple1(SignalDAQPeak,bins,title="Signal DAQ  peak in ADC (25 ns)",xlabel = "adc",
            save=CPAR['saveHistos'],filename='SignalDAQPeak.png', filepath=filepath)

        bins = hbins(SignalDAQOffMax, nsigma=5, nbins=20)
        HSimple1(SignalDAQOffMax,bins,title="Signal DAQ MAX ADC (25 ns)",xlabel = "adc",
            save=CPAR['saveHistos'],filename='SignalDAQOffMax.png', filepath=filepath)

        bins = hbins(SignalDAQOffMin, nsigma=5, nbins=20)
        HSimple1(SignalDAQOffMin,bins,title="Signal DAQ MIN ADC (25 ns)",xlabel = "adc",
            save=CPAR['saveHistos'],filename='SignalDAQOffMin.png', filepath=filepath)

        
        #plt.show()
def PlotSignals(DSGN,CPAR):
    """
    Plot signals
    """

    #spe
    if CPAR['plot_spe'] == True:
        if CPAR['plot_I']:
            print DSGN['spe_I']
            DSGN['spe_I'].PlotSignal(title= "SPE after PMT", units='muA')
        
        if CPAR['plot_V']:
            print DSGN['spe_V']
            DSGN['spe_V'].PlotSignal(title= "SPE (V) after PMT ", units='mV') 
                                      
        
        if CPAR['plot_ADC']:
            print DSGN['spe_ADC']
            DSGN['spe_ADC'].PlotSignal(title= "SPE (ADC) after PMT", units='adc') 
                                       
    
    if CPAR['plot_spe_fee']:


        print DSGN['spe_FEE']
        DSGN['spe_FEE_NN'].PlotSignal(title= "SPE after FEE NN" , units='mV') 
                                      
        DSGN['spe_FEE'].PlotSignal(title= "SPE after FEE" , units='mV')
                                    

        print DSGN['spe_DAQ']
        DSGN['spe_DAQ'].PlotSignal(title= "SPE after DAQ" , units='adc')


    # s1 and s2 MC average of 12 channels
    if CPAR['plot_s1_mc'] == True:
        print DSGN['s1_MC']
        DSGN['s1_MC'].PlotSignal(title= "S1 MC: positive downscaled") 
        

    if CPAR['plot_s2_mc'] == True:
        print DSGN['s2_MC']
        DSGN['s2_MC'].PlotSignal(title= "S2 MC: positive downscaled") 
       

    #s1 and s2 signal from PMT
    if CPAR['plot_s1_pmt'] == True:
        if CPAR['plot_I']:
            print DSGN['s1_PMT']
            DSGN['s1_PMT'].PlotSignal(title= "S1 after PMT: pos downscaled", units='muA') 
                                         
        
        if CPAR['plot_V']:
            print DSGN['s1_PMT_V']
            DSGN['s1_PMT_V'].PlotSignal(units='mV', 
                                        title= "S1 (V) after PMT: pos downscaled" )
        
        if CPAR['plot_ADC']:
            print DSGN['s1_PMT_ADC']
            
            DSGN['s1_PMT_ADC'].PlotSignal(units='adc', 
                                          title= "S1 (ADC) after PMT: pos downscaled")

    if CPAR['plot_s2_pmt'] == True:
        if CPAR['plot_I']:
            print DSGN['s2_PMT']
            
            DSGN['s2_PMT'].PlotSignal(units='muA', 
                                      title= "S2 after PMT: pos downscaled")
       
        if CPAR['plot_V']:
            print DSGN['s2_PMT_V']
            
            DSGN['s2_PMT_V'].PlotSignal(units='mV',
                                        title= "S2 (V) after PMT: pos downscaled")
        
        if CPAR['plot_ADC']:
            print DSGN['s2_PMT_ADC']
        
            DSGN['s2_PMT_ADC'].PlotSignal(units='adc',
                                          title= "S2 (ADC) after PMT: pos downscaled")
        
    if CPAR['plot_signal_pmt'] == True:

        print DSGN['signal_PMT_ADC']
        #
        DSGN['signal_PMT_ADC'].PlotSignal(units='adc',
                                          title= "Signal (ADC) after PMT: pos downscaled")
        
    #s1 and s2 signal from FEE
    if CPAR['plot_s1_fee'] == True:
        print DSGN['s1_FEE_NN']
        DSGN['s1_FEE_NN'].PlotSignal(units='mV', title= "S1 after FEE NN" )
        
        print DSGN['s1_FEE']
        DSGN['s1_FEE'].PlotSignal(units='mV', title= "S1 after FEE " )
        
        print DSGN['s1_DAQ']
        DSGN['s1_DAQ'].PlotSignal(units='adc', title= "S1 after DAQ " )

    if CPAR['plot_s2_fee'] == True:
        print DSGN['s2_FEE']
        DSGN['s2_FEE'].PlotSignal(units='mV', title= "S2 after FEE" )
        
        print DSGN['s2_DAQ']
        DSGN['s2_DAQ'].PlotSignal(units='adc', title= "S2 after DAQ" )

    if CPAR['plot_signal_fee'] == True:
        print DSGN['signal_DAQ']
        DSGN['signal_DAQ'].PlotSignal(units='adc', title= "Signal after DAQ" )
        print DSGN['signal_DAQ_off']
        DSGN['signal_DAQ_off'].PlotSignal(units='adc', 
                                            title= "Signal after DAQ offset" )
        
    #recovered
    if CPAR['plot_R'] == True:
        print DSGN['signal_R']
        DSGN['signal_R'].PlotSignal(units='adc', title= "Signal Recovered" )
        
        # DSGN['s2_R'].PlotSignal(positive = False, title= "S2 Recovered" )
        # print DSGN['s2_R']





