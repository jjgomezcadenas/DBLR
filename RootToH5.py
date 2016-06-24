"""
This module converts a root file containing histograms in an h5 file containing arrays
"""

from ROOT import TFile,gDirectory, gROOT
import RootUtil as RU 
import numpy as np
import tables as tb
import sys
from Util import wait


histo_bins = 1200000
nhistos = 100
path = "/Users/jjgomezcadenas/Documents/Development/NEXT/data/1ns/"

if __name__ == '__main__':

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

    dataset = ["bb0nu.1ns.15bar"]


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
        
        h5file = tb.open_file(path+h5file_name[ifile], mode = "w", title = h5file_name[ifile])
        root = h5file.root

        file = TFile( path+root_file_name[ifile] )
        keys=file.GetListOfKeys()
        rdir = keys[0]
        rdir.Print()

        wait()

        file.cd(rdir.GetName())
        khistos =gDirectory.GetListOfKeys()
        cnt = np.zeros(histo_bins, dtype=np.double)

        atom = tb.Atom.from_dtype(np.dtype('double'))
        filters = tb.Filters(complevel=5, complib='zlib')

        for khisto in khistos:
            histo =  gROOT.FindObject(khisto.GetName())
            name = histo.GetName()
            nbins= histo.GetNbinsX() 

            print "reading histogram =%s, with %d bins"%(name,nbins)
            if nbins != histo_bins:
                print "error: read %d bins, expect %d bins"%(nbins,histo_bins)
                sys.exit()

            RU.GetHistoContent(cnt, nbins, histo)
            hcnt = h5file.create_carray(root, name, atom, cnt.shape, filters=filters)
            hcnt[:] = cnt

        print(h5file)
        h5file.close()

    



    