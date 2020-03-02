"""
Makes an array (masked or not) from a filterbank file. Also applies bandpass correction using an offpulse region selected using RFI_offpulse.py.
Kenzie Nimmo 2020
"""

import filterbank
import numpy as np
import pickle
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt

def filterbank_to_np(filename, maskfile=None, dm=None, bandpass=False, offpulse=None, nbins=6):
    fil = filterbank.filterbank(filename)
    total_N = fil.number_of_samples
    spec=fil.get_spectra(0,total_N)
    if dm!=None:
        spec.dedisperse(dm, padval='mean')
    arr = np.array([spec[i] for i in xrange(fil.header['nchans'])])
    t_samp = fil.header['tsamp']
    if maskfile!=None:
        amaskfile = pickle.load(open(maskfile,'rb'))
        amask=[int(i) for i in amaskfile]
        vmin = np.amin(arr)
        arr[amask,:]=vmin-100
        mask = arr<vmin-50
        arr = np.ma.masked_where(mask==True,arr)
    arr=np.flip(arr,0)
    if bandpass==True and offpulse!=None:
        arr=bp(filename,maskfile,nbins,offpulse)
    return arr


def bp(filename,maskfile,nbins,offpulsefile,kind='slinear'):
    """
    Uses off pulse data (identified using interactive RFI_masker.py) to normalise the spectrum and
    calibrate the bandpass

    Note the spline fit is not working. It introduces weird artefacts. Currently just simply dividing without smoothing.
    """
    fil = filterbank.filterbank(filename)
    total_N = fil.number_of_samples
    nchans = fil.header['nchans']
    spec=fil.get_spectra(0,total_N)
    arr = np.array([spec[i] for i in xrange(fil.header['nchans'])])
    offpulse=pickle.load(open(offpulsefile,'rb'))
    spec = np.mean(arr[:,offpulse],axis=1)
    means=np.mean(spec[1:-1].reshape(-1, nbins), axis=1)
    meanint = np.append(spec[1],means)
    meanspec = np.append(meanint,spec[-1])
    x1 = np.append(0,np.arange(1,nchans-2,nbins))
    x = np.append(x1,nchans)
    cs = interp(x, meanspec, kind=kind)
    fit = cs(np.arange(0,nchans,1))
    if maskfile!=None:
        amaskfile = pickle.load(open(maskfile,'rb'))
        amask=[int(i) for i in amaskfile]
        amask = np.array(amask)
        vmin=np.amin(fit)
        fit[amask]=vmin-100
        mask = fit<vmin-50
        maskfit = np.ma.masked_where(mask==True,fit)
        vmin = np.amin(arr)
        arr[amask,:]=vmin-100
        maskar = arr<vmin-50
        arr = np.ma.masked_where(maskar==True,arr)
    else:
        maskfit = fit
    arr2=arr.copy()
    for i in range(arr.shape[1]):
        arr2[:,i]/=spec
        arr[:,i]/=maskfit

    #plt.plot(np.mean(arr,axis=1),'b')
    #plt.plot(np.mean(arr2,axis=1),'r')
    #plt.show()

    arr2-=np.mean(arr2)
    return arr2
