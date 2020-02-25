import filterbank
import numpy as np
import rfifind


def filterbank_to_np(filename, maskfile, total_N):
    fil = filterbank.filterbank(filename)
    spec=fil.get_spectra(0,total_N)
    arr = np.array([spec[i] for i in xrange(fil.header['nchans'])])
    t_samp = fil.header['tsamp']
    if maskfile:
        mask = rfifind.rfifind(maskfile)
        maskchans = mask.mask_zap_chans
        maskchans = fil.header['nchans']-1-np.array(list(maskchans))
        arr[maskchans,:]=0
    arr=np.flip(arr,0)
    return arr

def burst_spectra(numpy_arr,burst_cent,burst_width):
    burst_begin = int(np.ceil(burst_cent-(burst_width/2.)))
    burst_end = int(np.ceil(burst_cent+(burst_width/2.)))
    spec = np.mean(numpy_arr[:,burst_begin:burst_end],axis=1)
    return spec

def bandpass(arr,t_res):
    """
    Uses off pulse data (excluding 20ms around the burst) to normalise the spectrum and
    calibrate the bandpass
    """
    offset = int(10e-3/t_res)
    spectrum = np.mean(arr,axis=1)
    good_bins = spectrum!=0
    profile = np.mean(arr,axis=0)
    off_times = np.ones(arr.shape[1],np.bool)
    bin_peak = np.where(profile==np.amax(profile))[0][0]
    off_times[np.int(bin_peak-offset):np.int(bin_peak+offset)]=False
    offspec = np.mean(arr[:,off_times],axis=1)
    for i in range(arr.shape[1]):
        arr[:,i]/=offspec
    arr[~good_bins,:]=0
    return arr
