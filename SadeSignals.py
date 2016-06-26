
"""
Module with classes and utilities for signal manipulation in SADE
"""
from CParam import *
from PlotUtil import *
import FEE as FE
import FEParam as FP


class SignalDef:
    """
    Defines the signal
    """   
    def __init__(self, name, s1_l, s1_r, s2_l, s2_r):
        """
        name is the name of the data set
        s1w =[s1_l,s1_r]: bounds of s1 signal
        s2w =[s2_l,s2_r]: bounds of s2 signal

        pulse vector has the same dimensins of time: units PE, muA, mV, adc counts
        Positive signal is searched for above threshold
        pp is the positive pulse
        
        """
        self.name = name
        self.h5file = name + '.h5'
        self.s1_l = s1_l
        self.s1_r = s1_r
        self.s2_l = s2_l
        self.s2_r = s2_r

        self.sgn_l = s1_l
        self.sgn_r = s2_r

        self.len_s2 = s2_r - s2_l
        self.len_s1 = s1_r - s1_l

        self.len_signal = self.sgn_r - self.sgn_l

    def __str__(self): 

        s= """
            ----name = %s -----
            h5 file = %s
            s1_l = %d, s1_r = %d
            s2_l = %d, s2_r = %d
            sgn_l = %d, sgn_r = %d
            len_s1 = %d
            len_s2 = %d
            len_signal = %d
             
            
        """%(self.name, self.hfile, 
            self.s1_l,self.s1_r,self.s2_l,self.s2_r,self.sgn_l,self.sgn_r,
            self.len_s1,self.len_s2,self.len_signal)
        return s


class Signal:
    """
    Characterizes the signals
    """   
    def __init__(self, name, time, pulse, stype='NDAC', threshold = 0, downscale=False):
        """
        defines a signal 
        stype is the signal type: 
        NDAC = No-DAQ signal (a signal not previously decimated, in bins of 1 ns)
        DAX = DAQ signal (a signal previously decimated in bins of 25 ns)
        
        """
        if len(time) != len(pulse):
            print """error: len(time) must be equal to len(pulse)
                     found: len(time) = %d, len(pulse) = %d
            """%(len(time),len(pulse))
            sys.exit()
        
        self.name = name
        self.stype = stype
        self.downscale = downscale
        
        if self.downscale and self.stype =='NDAC':
            self.td, self.pd = FE.DownScaleSignal(time, pulse, int(FP.time_bin))
        else:
            self.td = time
            self.pd = pulse

        #zero-supression
        if threshold > 0:
            self.tp, self.pp  = DB.FindSignalAboveThr(self.td, self.pd, threshold = threshold)
        else:
            self.tp = self.td
            self.pp = self.pd

        self.time = self.tp
        self.pulse = self.pp

        if self.downscale or self.stype !='NDAC':
            self.tl = len(self.time)*FP.time_DAQ
            self.area = np.sum(self.pulse)*FP.time_DAQ
        else:
            self.tl = len(self.time)
            self.area = np.sum(self.pulse)


    def Length(self, units='ns'):
        """
        Return the length of the signal in units
        """
        xu = 1.*eval(units)
        return self.tl/xu

    def Area(self, units='pes'):
        """
        Return the area of the signal
        """
        xu = 1.*eval(units)
        return self.area/xu
        
    def Peak(self, units='pes'):
        """
        Return the peak of the signal (maximum value)
        """
        xu = 1.*eval(units)
        return np.amax(self.pulse)/xu

    def Minimum(self, units='pes'):
        """
        Return the minimum of the signal 
        """
        xu = 1.*eval(units)
        return np.amin(self.pulse)/xu

    def Maximum(self, units='pes'):
        """
        Return the maximum (Peak) of the signal 
        """
        
        return self.Peak(units)

    def Avg(self, positive=False, units='pes'):
        """
        Return the average height of the signal
        """
        xu = 1.*eval(units)
        return np.average(self.pulse)/xu

    def PlotSignal(self, title = 'signal',tunits = 'ns', units='pes'):
        """
        Plot signal 
        """
        xu = 1.*eval(units)
        tu = 1.*eval(tunits)
        plot_signal(self.time/tu,self.pulse/xu, 
                title = title, 
                signal_start=self.time[0]/tu, signal_end=self.time[-1]/tu, units=units)

    def __str__(self):    
        s= """
            Signal %s:

            length = %7.2f ns 
            area = %7.2f 
            peak (maximum) = %7.2f 
            minimum = %7.2f 
            avg = %7.2f 
            
        """%(self.name, self.Length(), self.Area(), 
            self.Peak(), self.Minimum(), self.Avg()
            )
        return s


def SignalAnalysis(CSGN,CPLOT, saveHistos=False, filepath ="./"):
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


    if CPLOT['Histograms']:

        
        bins = hbins(s1Area, nsigma=5, nbins=20)
        HSimple1(s1Area,bins,title="S1 area in PES (1ns)",xlabel = "pes",
            save=saveHistos,filename='s1Area.png', filepath=filepath)

        bins = hbins(s2Area, nsigma=5, nbins=20)
        HSimple1(s2Area,bins,title="S2 area (PES)",xlabel = "pes",
            save=saveHistos,filename='s2Area.png', filepath=filepath)

        bins = hbins(s2Length, nsigma=5, nbins=20)
        HSimple1(s2Length,bins,title="S2 length (mus)",xlabel = "mus",
            save=saveHistos,filename='s2Length.png', filepath=filepath)

        bins = hbins(s2Peak, nsigma=5, nbins=20)
        HSimple1(s2Peak,bins,title="S2 peak in PES (25 ns)",xlabel = "pes",
            save=saveHistos,filename='s2Peak.png', filepath=filepath)

        bins = hbins(s2Avg, nsigma=5, nbins=20)
        HSimple1(s2Avg,bins,title="S2 avg in PES (25 ns)",xlabel = "pes",
            save=saveHistos,filename='s2Avg.png', filepath=filepath)

        bins = hbins(SignalRecPeak, nsigma=5, nbins=20)
        HSimple1(SignalRecPeak,bins,title="Signal R  peak in ADC (25 ns)",xlabel = "adc",
            save=saveHistos,filename='SignalRecPeak.png', filepath=filepath)
        
        bins = hbins(SignalDAQPeak, nsigma=5, nbins=20)
        HSimple1(SignalDAQPeak,bins,title="Signal DAQ  peak in ADC (25 ns)",xlabel = "adc",
            save=saveHistos,filename='SignalDAQPeak.png', filepath=filepath)

        bins = hbins(SignalDAQOffMax, nsigma=5, nbins=20)
        HSimple1(SignalDAQOffMax,bins,title="Signal DAQ MAX ADC (25 ns)",xlabel = "adc",
            save=saveHistos,filename='SignalDAQOffMax.png', filepath=filepath)

        bins = hbins(SignalDAQOffMin, nsigma=5, nbins=20)
        HSimple1(SignalDAQOffMin,bins,title="Signal DAQ MIN ADC (25 ns)",xlabel = "adc",
            save=saveHistos,filename='SignalDAQOffMin.png', filepath=filepath)

        
        #plt.show()
def PlotSignals(DSGN,CPLOT):
    """
    Plot signals
    """

    #spe
    if CPLOT['plot_spe'] == True:
        if CPLOT['plot_I']:
            print DSGN['spe_I']
            DSGN['spe_I'].PlotSignal(title= "SPE after PMT", units='muA')
        
        if CPLOT['plot_V']:
            print DSGN['spe_V']
            DSGN['spe_V'].PlotSignal(title= "SPE (V) after PMT ", units='mV') 
                                      
        
        if CPLOT['plot_ADC']:
            print DSGN['spe_ADC']
            DSGN['spe_ADC'].PlotSignal(title= "SPE (ADC) after PMT", units='adc') 
                                       
    
    if CPLOT['plot_spe_fee']:


        print DSGN['spe_FEE']
        DSGN['spe_FEE_NN'].PlotSignal(title= "SPE after FEE NN" , units='mV') 
                                      
        DSGN['spe_FEE'].PlotSignal(title= "SPE after FEE" , units='mV')
                                    

        print DSGN['spe_DAQ']
        DSGN['spe_DAQ'].PlotSignal(title= "SPE after DAQ" , units='adc')


    # s1 and s2 MC average of 12 channels
    if CPLOT['plot_s1_mc'] == True:
        print DSGN['s1_MC']
        DSGN['s1_MC'].PlotSignal(title= "S1 MC: positive downscaled") 
        

    if CPLOT['plot_s2_mc'] == True:
        print DSGN['s2_MC']
        DSGN['s2_MC'].PlotSignal(title= "S2 MC: positive downscaled") 
       

    #s1 and s2 signal from PMT
    if CPLOT['plot_s1_pmt'] == True:
        if CPLOT['plot_I']:
            print DSGN['s1_PMT']
            DSGN['s1_PMT'].PlotSignal(title= "S1 after PMT: pos downscaled", units='muA') 
                                         
        
        if CPLOT['plot_V']:
            print DSGN['s1_PMT_V']
            DSGN['s1_PMT_V'].PlotSignal(units='mV', 
                                        title= "S1 (V) after PMT: pos downscaled" )
        
        if CPLOT['plot_ADC']:
            print DSGN['s1_PMT_ADC']
            
            DSGN['s1_PMT_ADC'].PlotSignal(units='adc', 
                                          title= "S1 (ADC) after PMT: pos downscaled")

    if CPLOT['plot_s2_pmt'] == True:
        if CPLOT['plot_I']:
            print DSGN['s2_PMT']
            
            DSGN['s2_PMT'].PlotSignal(units='muA', 
                                      title= "S2 after PMT: pos downscaled")
       
        if CPLOT['plot_V']:
            print DSGN['s2_PMT_V']
            
            DSGN['s2_PMT_V'].PlotSignal(units='mV',
                                        title= "S2 (V) after PMT: pos downscaled")
        
        if CPLOT['plot_ADC']:
            print DSGN['s2_PMT_ADC']
        
            DSGN['s2_PMT_ADC'].PlotSignal(units='adc',
                                          title= "S2 (ADC) after PMT: pos downscaled")
        
    if CPLOT['plot_signal_pmt'] == True:

        print DSGN['signal_PMT_ADC']
        #
        DSGN['signal_PMT_ADC'].PlotSignal(units='adc',
                                          title= "Signal (ADC) after PMT: pos downscaled")
        
    #s1 and s2 signal from FEE
    if CPLOT['plot_s1_fee'] == True:
        print DSGN['s1_FEE_NN']
        DSGN['s1_FEE_NN'].PlotSignal(units='mV', title= "S1 after FEE NN" )
        
        print DSGN['s1_FEE']
        DSGN['s1_FEE'].PlotSignal(units='mV', title= "S1 after FEE " )
        
        print DSGN['s1_DAQ']
        DSGN['s1_DAQ'].PlotSignal(units='adc', title= "S1 after DAQ " )

    if CPLOT['plot_s2_fee'] == True:
        print DSGN['s2_FEE']
        DSGN['s2_FEE'].PlotSignal(units='mV', title= "S2 after FEE" )
        
        print DSGN['s2_DAQ']
        DSGN['s2_DAQ'].PlotSignal(units='adc', title= "S2 after DAQ" )

    if CPLOT['plot_signal_fee'] == True:
        print DSGN['signal_DAQ']
        DSGN['signal_DAQ'].PlotSignal(units='adc', title= "Signal after DAQ" )
        print DSGN['signal_DAQ_off']
        DSGN['signal_DAQ_off'].PlotSignal(units='adc', 
                                            title= "Signal after DAQ offset" )
        
    #recovered
    if CPLOT['plot_R'] == True:
        print DSGN['signal_R']
        DSGN['signal_R'].PlotSignal(units='adc', title= "Signal Recovered" )
        
        # DSGN['s2_R'].PlotSignal(positive = False, title= "S2 Recovered" )
        # print DSGN['s2_R']



