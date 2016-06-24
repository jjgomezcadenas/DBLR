"""
This module handles ROOT IO
"""

from ROOT import TFile,gROOT
from ROOT import gDirectory
import numpy as np
import tables as tb
from PlotUtil import *
from Util import *


class HNP1:
    """
    stores a root histogram into numpy arrays
    """
    def __init__(self, histo):
        self.histo = histo
        self.name = histo.GetName()
        self.nbins= histo.GetNbinsX() 
        self.cnt = np.zeros(self.nbins, 'Float32')
        for i in xrange(0, self.nbins):
            self.cnt[i] = histo.GetBinContent (i)

    def __str__(self):
        
        s= """
        hnp1:
        name = %s, nbins = %d,   
        """%(self.name, self.nbins)
        return s

    
def get_histos(file_path, number_of_histos):
    HISTOS=[]

    file = TFile( file_path )
    keys=file.GetListOfKeys()
    rdir = keys[0]
    rdir.Print()
    file.cd(rdir.GetName())
    khistos =gDirectory.GetListOfKeys()

    n = 0
    for khisto in khistos:
        print "storing histo %s in list"%khisto.GetName()
        histo =  gROOT.FindObject(khisto.GetName())
        hnp = HNP1(histo)
        HISTOS.append(hnp)
        n+=1
        if n > number_of_histos: break
    return HISTOS
        
if __name__ == '__main__':

    nhistos = 100

    path = "/Users/jjgomezcadenas/Documents/Development/NEXT/data/1ns/"

    # dataset = ["electrons.1275keV.1ns",                 
    #         "electrons.2615keV.1ns",                
    #         "electrons.511keV.1ns",               
    #         "Kr83.z10cm.1ns",                      
    #         "Kr83.z15cm.1ns",                       
    #         "Kr83.z25cm.1ns",                       
    #         "Kr83.z5cm.1ns",
    #         "bb0nu.1ns",
    #          "electrons.1275keV.1ns.15bar"
    #          "electrons.2615keV.1ns.15bar",
    #          "electrons.511keV.1ns.15bar"]

    dataset = ["electrons.2615keV.1ns.15bar"]


    root_file_name = []
    h5file_name = []

    for d in dataset:
        root_file_name.append(d+".histos.root")
        h5file_name.append(d+".h5")
    
    for ifile in xrange(0,len(dataset)):
    
        print """
        Converting root histograms to h5
        root file = %s
        h5 file = %s
        """%(root_file_name[ifile],h5file_name[ifile])

        wait()
        histoList = get_histos(path+root_file_name[ifile], nhistos)
        h5file = tb.open_file(path+h5file_name[ifile], mode = "w", title = h5file_name[ifile])
        root = h5file.root
            
        for h0 in histoList:
        
            atom = tb.Atom.from_dtype(np.dtype('Float32'))
            filters = tb.Filters(complevel=5, complib='zlib')
            #event = h5file.createGroup(root, h0.name)
        
            #hcbin = h5file.create_carray(event, 'cbin', atom, h0.cbin.shape, filters=filters)
            hcnt = h5file.create_carray(root, h0.name, atom, h0.cnt.shape, filters=filters)
            #hcbin[:] = h0.cbin
            hcnt[:] = h0.cnt

        print(h5file)
        h5file.close()

    



    