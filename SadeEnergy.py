"""
Energy classes
"""
from Util import *
import numpy as np

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
        self.renergy25ns = MeanWOL(self.epes25ns,nmean=self.nmean, nsigma=self.nsigma)
        self.renergy1ns = MeanWOL(self.epes1ns,nmean=self.nmean, nsigma=self.nsigma)
        self.tenergy25ns = MeanWOL(self.q25ns,nmean=self.nmean, nsigma=self.nsigma)
        self.tenergy1ns = MeanWOL(self.q1ns,nmean=self.nmean, nsigma=self.nsigma)
        self.RTtoT25ns = MeanWOL(self.epes25nsR,nmean=self.nmean, nsigma=self.nsigma)
        self.RTtoT1ns = MeanWOL(self.epes1nsR,nmean=self.nmean, nsigma=self.nsigma)

    def AdcToPES(self, ebin='25ns'):
        """
        Returns the ratio ADC to PES
        """
        if ebin == '25ns':
            return np.average(self.adcToPes25ns)
        else:
            return np.average(self.adcToPes1ns)

    def RecoveredEnergyInPES(self, ebin='25ns'):
        """
        Returns the recovered (reconstructed) energy in a window
        of +- nsigma around the mean computed with the nmean first values
        """
        if ebin == '25ns':
            return self.renergy25ns
        else:
            return self.renergy1ns

    def RecoveredEnergyEff(self, ebin='25ns'):
        """
        Returns the fraction of events with reconstructed energy in a window
        of +- nsigma around the mean computed with the nmean first values
        """
        if ebin == '25ns':
            return len(self.renergy25ns)*1./len(self.epes25ns)
        else:
            return len(self.renergy1ns)*1./len(self.epes1ns)

    def RatioTrueToRecovered(self, ebin='25ns'):
        """
        Compute the difference ratio 2*abs(x1-x2)/(x1+x2)
        between true and reconstructed energy in a window
        of +- nsigma around the mean computed with the nmean first values
        """
        if ebin == '25ns':
            return  self.RTtoT25ns
        else:
            return self.RTtoT1ns

    def RatioTrueToRecoveredEff(self, ebin='25ns'):
        """
        Returns the fraction of events with diff ratio in a window
        of +- nsigma around the mean computed with the nmean first values
        """
        if ebin == '25ns':
            return len(self.RTtoT25ns)*1./len(self.epes25nsR)
        else:
            return len(self.RTtoT1ns)*1./len(self.epes1nsR)

    def TrueEnergyInPES(self, ebin='25ns'):
        """
        Returns the true energy in a window
        of +- nsigma around the mean computed with the nmean first values
        """
        if ebin == '25ns':
            return self.tenergy25ns
        else:
            return self.tenergy1ns

    def TrueEnergyEff(self, ebin='25ns'):
        """
        Returns the fraction of events with true energy in a window
        of +- nsigma around the mean computed with the nmean first values
        """
        if ebin == '25ns':
            return len(self.tenergy25ns)*1./len(self.q25ns)
        else:
            return len(self.tenergy1ns)*1./len(self.q1ns)
            

    def __str__(self):    
        s= """
           Energy
            adcToPES = %7.2f

            True Energy (1 ns) in pes
            avg = %7.2f, std = %7.2f, 

            True Energy in Window (1 ns) in pes
            eff = %7.2f
            avg = %7.2f, std = %7.2f, sigma_E/E (FWHM) = %7.2f

            True Energy (25 ns) in pes
            avg = %7.2f, std = %7.2f

            True Energy in Window (25 ns) in pes
            eff = %7.2f
            avg = %7.2f, std = %7.2f, sigma_E/E (FWHM) = %7.2f

            Recovered Energy in Window (25 ns) in pes
            eff = %7.2f
            avg = %7.2f, std = %7.2f, sigma_E/E (FWHM) = %7.2f

            Ratio True to Recovered Energy in Window (25 ns) in pes
            eff = %7.2f
            avg = %7.2g, std = %7.2g 

           
            """%(self.AdcToPES(ebin='25ns'),
                np.average(self.q1ns),np.std(self.q1ns),
                self.TrueEnergyEff(ebin='1ns'),
                np.average(self.TrueEnergyInPES(ebin='1ns')),
                np.std(self.TrueEnergyInPES(ebin='1ns')),
                ResFWHM(np.std(self.TrueEnergyInPES(ebin='1ns')),
                        np.average(self.TrueEnergyInPES(ebin='1ns'))),
                np.average(self.q25ns),np.std(self.q25ns),
                self.TrueEnergyEff(ebin='25ns'),
                np.average(self.TrueEnergyInPES(ebin='25ns')),
                np.std(self.TrueEnergyInPES(ebin='25ns')),
                ResFWHM(np.std(self.TrueEnergyInPES(ebin='25ns')),
                        np.average(self.TrueEnergyInPES(ebin='25ns'))),
                self.RecoveredEnergyEff(ebin='25ns'),
                np.average(self.RecoveredEnergyInPES(ebin='25ns')),
                np.std(self.RecoveredEnergyInPES(ebin='25ns')),
                ResFWHM(np.std(self.RecoveredEnergyInPES(ebin='25ns')),
                        np.average(self.RecoveredEnergyInPES(ebin='25ns'))),
                self.RatioTrueToRecoveredEff(ebin='25ns'),
                np.average(self.RatioTrueToRecovered(ebin='25ns')),
                np.std(self.RatioTrueToRecovered(ebin='25ns')))
        return s



    