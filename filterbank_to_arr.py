"""
Makes an array (masked or not) from a filterbank file. Also applies bandpass correction using an offpulse region selected using RFI_offpulse.py.
"""

import filterbank
import numpy as np
import pickle
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt

def filterbank_to_np(filename, dm=None, maskfile=None, bandpass=False, offpulse=None, nbins=6):
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

def fits_to_np(filename, dm=None, maskfile=None, bandpass=False, offpulse=None, n\
bins=6,AO=False):
    fits=psrfits.PsrfitsFile(filename)
    total_N=fits.specinfo.N
    t_samp=fits.specinfo.dt
    if AO==True:
        peak_bin=(total_N/10.)*2 #the chopped filterbank is 10 subints (2 before \
the burst and 8 after)
        #+-0.1seconds
        begin_bin=int(peak_bin-(0.1/t_samp)) #number of bins
        end_bin=int(peak_bin+(0.1/t_samp))

        spec=fits.get_spectra(begin_bin,end_bin)
    else: spec=fits.get_spectra(0,total_N)
    if dm!=None:
        spec.dedisperse(dm, padval='mean')
    arr = np.array([spec[i] for i in xrange(fits.specinfo.num_channels)])

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

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y



def bp(filename,maskfile,nbins,offpulsefile,kind='slinear'):
    """
    Uses off pulse data (identified using interactive RFI_masker.py) to normalise the spectrum and
    calibrate the bandpass
    """
    fil = filterbank.filterbank(filename)
    total_N = fil.number_of_samples
    nchans = fil.header['nchans']
    spec=fil.get_spectra(0,total_N)
    arr = np.array([spec[i] for i in xrange(fil.header['nchans'])])
    offpulse=pickle.load(open(offpulsefile,'rb'))
    spec = np.mean(arr[:,offpulse],axis=1)
    # smooth the bandpass spectrum
    spec_smooth = smooth(spec,window_len=10)

    plt.plot(spec)
    plt.plot(spec_smooth)
    plt.show()

    if maskfile!=None:
        amaskfile = pickle.load(open(maskfile,'rb'))
        amask=[int(i) for i in amaskfile]
        maskarr=np.zeros_like(arr)
        maskarr[amask]=1
        arr = np.ma.masked_where(maskarr==True,arr)
    arr2=arr.copy()
    for i in range(arr.shape[1]):
        arr2[:,i]/=spec_smooth
        #arr[:,i]/=maskfit

    arr2-=np.mean(arr2)
    plt.plot(np.mean(arr2,axis=1))
    plt.show()
    return arr2
