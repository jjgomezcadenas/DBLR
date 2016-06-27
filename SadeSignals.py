
"""
Module with classes and utilities for signal manipulation in SADE
"""
from CParam import cython_dblr

if cython_dblr == False:
    import DBLR as DB
else:
    import PxDBLR as DB
    
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
            self.tp, self.pp  = DB.FindSignalAboveThr(self.td, self.pd, 
                                                      threshold = threshold)
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


