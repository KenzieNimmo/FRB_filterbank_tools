"""
Kenzie Nimmo 2020
"""

import numpy as np
import filterbank_to_arr
from fit_repeater_bursts import fit_my_smudge
import filterbank
import os
import pickle
import matplotlib.pyplot as plt
import sys 
import re 

if __name__ == '__main__':
    path = './'
    IDs_ordered = ['4']
    pklfilename = sys.argv[1]
    """
    pklfilename['array_corrected']
    pklfilename['array_uncorrected']
    pklfilename['mask']
    pklfilename['t_samp']
    pklfilename['freqs']
    """
    bursts = pickle.load(open(pklfilename,'rb'))
    picklename = re.search('(.*).pkl',pklfilename).group(1)
    
    for burst in IDs_ordered:
        burst=str(burst)
        arr = bursts['array_corrected']
        #arr=filterbank_to_arr.bp('burst4.fil','burst4_mask.pkl',6,'burst4_offpulse_time.pkl') #temporary step that will be accounted for in RFI_offpulse.py
        prof = np.mean(arr,axis=0)
        spec = np.mean(arr,axis=1)
        total_N=len(prof)
        t_samp = bursts['t_samp']

        #fil=filterbank.filterbank('burst4.fil')
        #properties of filterbank file
        #freqs = np.flip(fil.frequencies)
        freqs = np.flip(bursts['freqs'])
        #f_res = (freqs[-1]-freqs[0])/(fil.header['nchans']-1)
 
        stimes = np.linspace(0,arr.shape[1],arr.shape[1])
        stimes*=t_samp
        #guess=[dmax, xmax, ymax, xwid, ywid, 0]
        if burst == "4":
            guess = [50,450,4000,4,8000,0]
        else: guess = []

        bin_times,fit, bin_center = fit_my_smudge(arr, stimes, freqs, guess=guess, doplot=True, basename=picklename)
        fits={}
        burst_properties={}
        fits['array_corrected']=bursts['array_corrected']
        fits['array_uncorrected']=bursts['array_uncorrected']                                                                            
        fits['mask']=bursts['mask']
        fits['t_samp']=bursts['t_samp']
        fits['freqs']=bursts['freqs']
        fits['centre_bin']=bin_times[0]
        fits['width_bin']=bin_times[1]
        
        burst_properties[burst] = fits 
    f=open("%s.pkl"%picklename, "wb")
    pickle.dump(burst_properties,f)
    f.close()
        
