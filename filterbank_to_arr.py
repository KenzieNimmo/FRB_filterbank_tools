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
    return arr, t_samp

def spectra(numpy_arr,burst_cent,burst_width):
    burst_begin = int(np.ceil(burst_cent-(burst_width/2.)))
    burst_end = int(np.ceil(burst_cent+(burst_width/2.)))
    spec = np.mean(numpy_arr[:,burst_begin:burst_end],axis=1)
    return spec
