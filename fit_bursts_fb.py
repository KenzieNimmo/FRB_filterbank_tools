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

if __name__ == '__main__':
    path = './'
    IDs_ordered = ["4"]
    bursts = {"4": {'fb': 'test.fil', 'mask': None}

    }

    for burst in IDs_ordered:
        filname=bursts[burst]['fb']
        maskname=bursts[burst]['mask']
        filfile=os.path.join(path,filname)
        maskfile=maskname
        if bursts[burst]['mask']!=None:
            maskfile=os.path.join(path,maskname)

        fil = filterbank.filterbank(filname)
        total_N=1000
        t_samp = fil.header['tsamp']

        arr = filterbank_to_arr.filterbank_to_np(filfile,maskfile,total_N)
        amask_file=None
        if amask_file!=None:
            amaskfile = np.loadtxt(amask_file)
            amask=[int(i) for i in amaskfile]
            arr[amask,:]=0
        arr_burst = filterbank_to_arr.bandpass(arr, t_samp) #bandpass correct
        
        
        #properties of filterbank file
        freqs = np.flip(fil.frequencies)
        f_res = (freqs[-1]-freqs[0])/(fil.header['nchans']-1)
        stimes = np.linspace(0,arr_burst.shape[1],arr_burst.shape[1])
        stimes*=t_samp
        #guess=[dmax, xmax, ymax, xwid, ywid, 0]
        if burst == "4":
            guess = [50,450,4000,4,8000,0]
        else: guess = []

        bin_times,fit, bin_center = fit_my_smudge(arr_burst, stimes, freqs, guess=guess, doplot=True, basename=filname)
        fits={}
        burst_properties={}
        fits['centre_bin']=bin_times[0]
        fits['width_bin']=bin_times[1]
        
        burst_properties[burst] = fits 
    f=open("bursts_peak_width.pkl", "wb")
    pickle.dump(burst_properties,f)
    f.close()
        
