"""
Kenzie Nimmo 2020
"""

import numpy as np
import filterbank_to_arr
from fit_repeater_bursts import fit_my_smudge
import filterbank
import os



if __name__ == '__main__':
    path = './'
    IDs_ordered = ["4"]
    bursts = {"4": {'fb': 'test.fil'}

    }

    for burst in IDs_ordered:
        filname=bursts[burst]['fb']
        filfile=os.path.join(path,filname)
        fil = filterbank.filterbank(filname)
        total_N=1000
        t_samp = fil.header['tsamp']

        arr = filterbank_to_arr.filterbank_to_np(filfile,None,total_N)
        amask_file="./mask.txt"
        if amask_file!=None:
            amaskfile = np.loadtxt(amask_file)
            amask=[int(i) for i in amaskfile]
            arr[amask,:]=0

        spec_bp_corr = filterbank_to_arr.bandpass(arr, t_samp) #bandpass correct

        #properties of filterbank file
        freqs = np.flip(fil.frequencies)
        f_res = (freqs[-1]-freqs[0])/(fil.header['nchans']-1)
        stimes = np.linspace(0,total_N,total_N)
        stimes*=t_samp
        #guess=[dmax, xmax, ymax, xwid, ywid, 0]
        if burst == "100":
            guess = [50,100,57,10,12,0]
        else: guess = []

        fit, bin_center = fit_my_smudge(arr, stimes, freqs, guess=guess, doplot=True, basename=filname)
        #fits[str(burst)] = fit[3:-1]
        #bin_centers = np.vstack([bin_centers,bin_center])
