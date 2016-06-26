"""
Energy calculations
"""
from PlotUtil import *
import numpy as np
from scipy import signal as SGN

class EnergyVectors:
    """
    Stores the data needed for energy calculation
    """   
    def __init__(self):
        """
        q1ns : S2 signal in pes (bins of 1 ns)
        q25ns : S2 signal in pes (bins of 25 ns) 
        adc1ns : S2 signal in adc (bins of 1 ns)
        adc25ns : S2 signal in adc (bins of 25 ns) 
        eadc: Reconstructed signal in adc (bins of 25 ns)
        """

        self.i = 0
        self.q1ns = []
        self.q25ns = []
        self.adc1ns = []
        self.adc25ns = []
        self.eadc = []
        self.epes1ns = [] 
        self.epes25ns = [] 
        self.qR = []
        self.adcToPes1ns = []
        self.adcToPes25ns = []
        self.adcToPesR = []
        self.epes1nsR = [] 
        self.epes25nsR = []

    def AddVectors(self,q_1ns,q_25ns,adc_1ns,adc_25ns,eadc):
        """
        Collects data
        """
        self.i+=1
        self.q1ns.append(q_1ns)
        self.q25ns.append(q_25ns)
        self.adc1ns.append(adc_1ns)
        self.adc25ns.append(adc_25ns)
        self.eadc.append(eadc)
        self.qR.append(differenceRatio(q_1ns,q_25ns))
        self.adcToPes1ns.append(adc_1ns/q_1ns)
        self.adcToPes25ns.append(adc_25ns/q_25ns)
        self.adcToPesR.append(differenceRatio(adc_1ns/q_1ns,adc_25ns/q_25ns))
        self.epes1ns.append(eadc/(adc_1ns/q_1ns))
        self.epes25ns.append(eadc/(adc_25ns/q_25ns))
        self.epes1nsR.append(differenceRatio(eadc/(adc_1ns/q_1ns),q_1ns))
        self.epes25nsR.append(differenceRatio(eadc/(adc_25ns/q_25ns),q_25ns))

    def __str__(self):    
        s= """
           Energy Vectors, event = %d

            Q: S2 MC (pes): 1 ns  = %7.2f
            Q: S2 MC (pes): 25 ns  = %7.2f
            Q: S2 PMT (adc): 1 ns  = %7.2f
            Q: S2 PMT (adc): 25 ns  = %7.2f
            QR: S2 MC 1/25 ns  = %7.2g
            energy (adc) = %7.2f
            adcToPes (1 ns) = %7.2f 
            adcToPes (25 ns) = %7.2f
            adcToPesR:  1/25 ns  = %7.2g
            energy (pes): 1 ns = %7.2f 
            energy (pes): 25 ns = %7.2f  
            Energy to S2 MC: difference ratio (1ns) = %7.2g
            Energy to S2 MC: difference ratio (25ns) = %7.2g
           
            """%(self.i,
                self.q1ns[self.i-1],
                self.q25ns[self.i-1],
                self.adc1ns[self.i-1],
                self.adc25ns[self.i-1],
                self.qR[self.i-1],
                self.eadc[self.i-1],
                self.adcToPes1ns[self.i-1],
                self.adcToPes25ns[self.i-1],
                self.adcToPesR[self.i-1],
                self.epes1ns[self.i-1],
                self.epes25ns[self.i-1],
                self.epes1nsR[self.i-1],
                self.epes25nsR[self.i-1])
        return s

class SadeEnergy:
    """
    Energy calculation
    """   
    def __init__(self, energyVectors, nmean=5, nsigma=5):

        self.q1ns = np.array(energyVectors.q1ns)
        self.q25ns = np.array(energyVectors.q25ns)
        self.adc1ns = np.array(energyVectors.adc1ns)
        self.adc25ns = np.array(energyVectors.adc25ns)
        self.eadc = np.array(energyVectors.eadc)
        self.qR = np.array(energyVectors.qR)
        self.adcToPes1ns = np.array(energyVectors.adcToPes1ns)
        self.adcToPes25ns = np.array(energyVectors.adcToPes25ns)
        self.adcToPesR = np.array(energyVectors.adcToPesR)
        self.epes1ns = np.array(energyVectors.epes1ns)
        self.epes25ns = np.array(energyVectors.epes25ns)
        self.epes1nsR = np.array(energyVectors.epes1nsR)
        self.epes25nsR = np.array(energyVectors.epes25nsR)
        self.nmean = nmean
        self.nsigma = nsigma


    def RecoveredEnergyInPES(self, ebin='25ns'):
        if ebin == '25ns':
            return MeanWOL(self.epes25ns,nmean=self.nmean, nsigma=self.nsigma)
        else:
            return MeanWOL(self.epes1ns,nmean=self.nmean, nsigma=self.nsigma)

    def RatioTrueToRecovered(self, ebin='25ns'):
        if ebin == '25ns':
            return MeanWOL(self.epes25nsR,nmean=self.nmean, nsigma=self.nsigma) 
        else:
            return MeanWOL(self.epes1nsR,nmean=self.nmean, nsigma=self.nsigma)

    def TrueEnergyInPES(self, ebin='25ns'):
        if ebin == '25ns':
            return self.q25ns
        else:
            return self.q1ns

    def __str__(self):    
        s= """
           Energy

            Q: S2 MC (pes): 1 ns: avg = %7.2f, std = %7.2f
            Q: S2 MC (pes): 25 ns: avg = %7.2f, std = %7.2f
            Q: S2 PMT (adc): 1 ns: avg = %7.2f, std = %7.2f
            Q: S2 PMT (adc): 25 ns: avg = %7.2f, std = %7.2f
            QR: S2 MC 1/25 ns: avg = %7.2g, std = %7.2g
            energy (adc) : avg = %7.2f, std = %7.2f
            adcToPes (1 ns) : avg = %7.2f, std = %7.2f
            adcToPes (25 ns) : avg = %7.2f, std = %7.2f
            adcToPesR:  1/25 ns  : avg = %7.2g, std = %7.2g
            energy (pes): 1 ns : avg = %7.2f, std = %7.2f
            energy (pes): 25 ns : avg = %7.2f, std = %7.2f  
            Energy to S2 MC: difference ratio (1ns) : avg = %7.2g, std = %7.2g
            Energy to S2 MC: difference ratio (25ns) : avg = %7.2g, std = %7.2g
           
            """%(np.average(self.q1ns),np.std(self.q1ns),
                np.average(self.q25ns),np.std(self.q25ns),
                np.average(self.adc1ns),np.std(self.adc1ns),
                np.average(self.adc25ns),np.std(self.adc25ns),
                np.average(self.qR),np.std(self.qR),
                np.average(self.eadc),np.std(self.eadc),
                np.average(self.adcToPes1ns),np.std(self.adcToPes1ns),
                np.average(self.adcToPes25ns),np.std(self.adcToPes25ns),
                np.average(self.adcToPesR),np.std(self.adcToPesR),
                np.average(self.epes1ns),np.std(self.epes1ns),
                np.average(self.epes25ns),np.std(self.epes25ns),
                np.average(self.epes1nsR),np.std(self.epes1nsR),
                np.average(self.epes25nsR),np.std(self.epes25nsR))
        return s


def EnergyAnalysis(SE,CPLOT, saveHistos=False, filepath ="./"):
    """
    Takes an instance of Sade Energy
    """

    print """
        Sade Energy = %s
    """%(SE)

    epes25ns = SE.RecoveredEnergyInPES(ebin='25ns')
    
    qr25ns = 1000*SE.RatioTrueToRecovered(ebin='25ns')
    qr1ns = 1000*SE.RatioTrueToRecovered(ebin='1ns')

    etrue25ns = SE.TrueEnergyInPES(ebin='25ns')
    etrue1ns = SE.TrueEnergyInPES(ebin='1ns')

    print """

        energy vector length = %d

        True Energy (MC only FANO): 25 ns
        mean = %7.2f std = %7.2f, sigma_E/E (FWHM) = %7.2f

        True Energy (MC only FANO): 1 ns
        mean = %7.2f std = %7.2f, sigma_E/E (FWHM) = %7.2f 

        Difference Ratio (x1000) = %7.2g

        Reconstructed energy (in PES, 25 ns includes electronics noise, recovery)
        mean = %7.2f std = %7.2f, sigma_E/E (FWHM) = %7.2f

        Difference Ratio True to Recovered (25 ns) x 1000 = %7.2f
        Difference Ratio True to Recovered (1 ns) x 1000 = %7.2f

    """%(len(epes25ns), 
         np.average(etrue25ns), np.std(etrue25ns),
         ResFWHM(np.std(etrue25ns),np.average(etrue25ns)),
         np.average(etrue1ns), np.std(etrue1ns),
         ResFWHM(np.std(etrue1ns),np.average(etrue1ns)),
         differenceRatio(np.average(etrue1ns),np.average(etrue25ns)),
         np.average(epes25ns), np.std(epes25ns),
         ResFWHM(np.std(epes25ns),np.average(epes25ns)),
         np.average(qr25ns),
         np.average(qr1ns)
         )

    if CPLOT['EnergyHistograms'] == True:

        bins = hbins(etrue25ns, nsigma=5, nbins=20)
        HSimple1(etrue25ns,bins,title="True energy in PES (25 ns)",xlabel = "pes",
            save=saveHistos,filename='etrue25ns.png', filepath=filepath)

        bins = hbins(epes25ns, nsigma=5, nbins=20)
        HSimple1(epes25ns,bins,title="Recovered energy in PES (25 ns)",xlabel = "pes",
            save=saveHistos,filename='epes25ns.png', filepath=filepath)

        bins = hbins(qr25ns, nsigma=5, nbins=20)
        HSimple1(qr25ns,bins,title="Rec ratio (25 ns)",xlabel = "",
            save=saveHistos,filename='qr25ns.png', filepath=filepath)


    