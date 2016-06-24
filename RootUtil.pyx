"""
A function to copy a root histogram into a numpy array
"""
def GetHistoContent(double[:] cnt, unsigned int nbins, histo):
	"""
	Copy the contents of a root histogram into array cnt
	"""
	cdef unsigned int i
	for i in range(nbins):
		cnt[i] = histo.GetBinContent(i)

	return cnt
