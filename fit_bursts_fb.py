"""
Kenzie Nimmo 2020
"""
import sys
sys.path.insert(1,'~/FRB_filterbank_tools')
import numpy as np
import filterbank_to_arr
from fit_repeater_bursts import fit_my_smudge
import filterbank
import os
import pickle
import matplotlib.pyplot as plt
import re 
import optparse
import barycenter_correct_fb

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [options] infile', \
                description="2D Gaussian fit to FRB data. Input the pickle file output from RFI_masker.py")
    parser.add_option('-g', '--guess', dest='guess', type='string', \
                      help="Guess for gaussian fit. time_peak:time_width:freq_peak:freq_width.", default=None)
    parser.add_option('-u', '--uncorr', dest='uncorr', action="store_true", \
                      help="If -u option is used, use the uncorrected array (otherwise use the masked+bandpass corrected array).", default=False)
    parser.add_option('-d', '--dm', dest='dm', type='float',help="Dispersion measure.", default=None)
    parser.add_option('-f', '--FRBname', dest='FRBname', type='string',help="FRB R name.", default=None)
    parser.add_option('-t', '--telescopename', dest='telescopename', type='string',help="Telescope used (Eff, CHIME, DSS43).", default=None)
    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        sys.exit(1)
    elif len(args)!=1:
        sys.stderr.write("Only one input file must be provided!\n")
    else:
        options.infile = args[-1]


    path = './'
    pklfilename = options.infile
    """
    pklfilename['array_corrected']
    pklfilename['array_uncorrected']
    pklfilename['mask']
    pklfilename['t_samp']
    pklfilename['freqs']
    """
    bursts = pickle.load(open(pklfilename,'rb'))
    picklename = re.search('(.*).pkl',pklfilename).group(1)
    
    
    arr = bursts['array_corrected']
    if options.uncorr==True:
        arr = bursts['array_uncorrected']
    #arr=filterbank_to_arr.bp('burst4.fil','burst4_mask.pkl',6,'burst4_offpulse_time.pkl') #temporary step that will be accounted for in RFI_offpulse.py
    prof = np.mean(arr,axis=0)
    spec = np.mean(arr,axis=1)
    total_N=len(prof)
    t_samp = bursts['t_samp']
    freqs = bursts['freqs']
    
    stimes = np.linspace(0,arr.shape[1],arr.shape[1])
    stimes*=t_samp
    #guess=[dmax, xmax, ymax, xwid, ywid, 0]
    if options.guess!=None:
        g=options.guess.split (':') #time_peak,time_width,freq_peak,freq_width
        g = [int(i) for i in g]
        guess = [50, g[0], g[2], g[1], g[3], 0]
    else: guess = []

    times,fit, bin_center = fit_my_smudge(arr, stimes, freqs, guess=guess, doplot=True, basename=picklename)
    #times is peak bin, width bin, peak time in seconds, width time in seconds.
    #if options.dm!=None:
        #nchan=len(spec)
        #chanwidth=(freqs[-1]-freqs[0]/nchan)
        #peak_time=barycenter_correct_fb.barycorr(bursts['tstart'],times[2],(freqs[-1]+chanwidth/2.),options.dm,FRB=str(options.FRBname),telescope=str(options.telescopename))
        #print(peak_time)


    fits={}

    fits['array_corrected']=bursts['array_corrected']
    fits['array_uncorrected']=bursts['array_uncorrected']                                                                            
    fits['mask']=bursts['mask']
    fits['t_samp']=bursts['t_samp']
    fits['freqs']=bursts['freqs']
    fits['centre_bin']=times[0]
    fits['width_bin']=np.abs(times[1])
        
    with open("%s.pkl"%picklename, 'wb') as f:    
        pickle.dump(fits,f)
    f.close()
        
