"""
KN 2020
"""
import os
import filterbank
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import rfifind
from scipy.optimize import curve_fit, leastsq
import filterbank_to_arr
import pickle
mm_to_in = 0.0393701

"""
Not working yet - first thing is to incorporate Laura's gaussian fit... Then next thing to try is to make an archive file of this data set
and use pazi (copy pazi channels here) and run through scint_bw_highfreq.py
That's the only differences I can think of so far.
"""

def shift(v, i, nchan):
        """
        function v by a shift i
        nchan is the number of frequency channels (to account for negative lag)
        """
        n = len(v)
        r = np.zeros(3*n)
        i+=nchan-1 #to account for negative lag
        i = int(i)
        r[i:i+n] = v
        return r

def autocorr(x, nchan, v=None):
    if v is None:
        v = np.ones_like(x)
    x = x.copy()
    x[v!=0] -= x[v!=0].mean()
    ACF = np.zeros_like(x)
    for i in range(len(x)):
        if i>0:
            m = shift(v,0,nchan)*shift(v,i,nchan)
            ACF[i-1] = np.sum(shift(x,0,nchan)*shift(x, i,nchan)*m)/np.sqrt(np.sum(shift(x, 0, nchan)**2*m)*np.sum(shift(x, i, nchan)**2*m))
    return ACF

def scint_bw(filename,maskfile,total_N,additional_masking=False):
    #For R3, we made the narrow peak one bin wide in time to achieve the best frequency resolution. But there would need to be an
    #additional step in here to fit a Gaussian to the burst profile to determine the peak and width.

    arr = filterbank_to_arr.filterbank_to_np(filename,maskfile,total_N)
    if additional_masking!=False:
        arr[amask,:]=0

    fil = filterbank.filterbank(filename)
    #inspect arr here using plt.imshow to see if you need to do additional masking
    prof = np.mean(arr,axis=0)
    bin_peak = np.where(prof==np.amax(prof))[0][0]

    #additional step to sum all the frequencies across the burst width required, in general. We only have one bin in time for R3.

    #get spectra for burst
    spec = filterbank_to_arr.burst_spectra(arr,bin_peak,6)

    #define an off-pulse region and subtract off pulse spectrum from burst spectrum
    t_res = fil.header['tsamp']
    offset = int(10e-3/t_res) #numbers of bins equivalent to 10ms.
    good_bins = spec!=0
    off_times = np.ones(arr.shape[1],np.bool)
    off_times[np.int(bin_peak-offset):np.int(bin_peak+offset)]=False
    offspec = np.mean(arr[:,off_times],axis=1)

    #normalize the spectrum
    spec_norm = spec/offspec
    spec_norm[~good_bins]=0 #making the zero spectrum values zero again after the normalization
    meanspec = np.mean(spec_norm)

    nchan = fil.header['nchans']
    channel_list = np.linspace(0,nchan-1,nchan)
    freqs = fil.frequencies
    f_res = (freqs[0]-freqs[-1])/nchan #MHz
    delta_f=np.linspace(0,nchan,nchan)*f_res

    good = np.ones(len(spec))
    good[~good_bins]=0 #good is an array of 1 and 0 where 1 is good and 0 is bad (so we ignore the bad)

    ACF=autocorr(spec_norm, nchan=nchan, v=good)
    plt.plot(ACF, 'r', label='ACF', drawstyle='steps-mid')

    return ACF, delta_f

def lorentzian_fit(ACF, delta_f, nbins_cent):
    """
    Fits a Lorentzian to the central nbins_cent bins of the ACF
    """
    # only focus on the central part of the ACF
    ACF_cent=ACF[1:nbins_cent]
    delta_cent = delta_f[1:nbins_cent]
    reversed_arr = ACF_cent[::-1]
    ACFc = np.append(reversed_arr,ACF_cent)
    maxACF=np.max(ACFc)
    reversed_arr = delta_cent[::-1]*-1
    delta_fc = np.append(reversed_arr, delta_cent)
    delta_fn = np.append(reversed_arr, np.array([-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05]))
    delta_fnew = np.append(delta_fn, delta_cent)

    #define the lorentzian function
    def lorentz(x,gamma, y0, c):
        return (y0)/(((x)**2)+gamma**2)+c

    uncert_region = ACF[192:1280]
    std = np.std(uncert_region)
    std_array = np.zeros_like(ACFc)
    std_array+=std

    popt, pcov = curve_fit(lorentz, delta_fc, ACFc, p0=[0.01,0.2,0], sigma=std_array)

    print "Scintillation BW:", popt[0], "MHz"

    return ACFc, delta_fc, lorentz(delta_fnew,*popt), delta_fnew, popt[0]

if __name__ == '__main__':

    #open the dictionary created in fit_bursts_fb.py
    with open('bursts_peak_width.pkl', 'rb') as f:
        bursts = pickle.load(f)

    path='./'
    IDs_ordered=['4']
    for burst in IDs_ordered:
        filname=bursts[burst]['fb']
        maskname=bursts[burst]['mask']
        filfile=os.path.join(path,filname)
        maskfile=maskname
        if bursts[burst]['mask']!=None:
            maskfile=os.path.join(path,maskname)

        #amask = np.linspace(6000, 7200, 1200)
        amask = np.loadtxt("./mask.txt")
        amask = [int(i) for i in amask]
        ACF, delta_f = scint_bw(filfile,maskfile,1000, additional_masking=amask)
        reversed_arr = delta_f[::-1]
        totdelta_f = np.append(-1*reversed_arr,delta_f[0:])
        ACFc, deltafc, lorentz, deltafnew,bw=lorentzian_fit(ACF, delta_f,60)

        fig = plt.figure(figsize=[183*mm_to_in,100*mm_to_in])
        grid = plt.GridSpec(2, 1, hspace=0.3, wspace=0.3)

        ax2 = fig.add_subplot(grid[0:1,:])
        ax2.plot(totdelta_f,np.append(ACF[::-1],ACF[0:]), color='darkorange', drawstyle='steps-mid')
        ax2.set_xticks([-128, -64, 0, 64, 128])
        ax2.set_yticks([0,0.1])
        ax2.set_yticklabels([r"$0$", r"$0.1$"])
        ax2.set_ylim([-0.02,0.11])

        ax4 = fig.add_subplot(grid[1:,:])
        ax4.plot(deltafc, ACFc, color='darkorange',drawstyle='steps-mid')
        ax4.plot(deltafnew, lorentz, color='darkgreen')
        ax4.axvline(x=bw,color='k',dashes=(5,2), lw=1.5 )
        ax4.set_xticks([-2,-1,0,1,2])
        ax4.set_yticks([0,0.06])
        ax4.set_yticklabels([r"$0$", r"$0.06$"])
        ax4.tick_params(axis='y', labelleft=True)
        ax4.set_xlim([-2,2])
        ax4.set_ylim([-0.035,0.073])

        fig.text(0.5, 0.03, r'Frequency Lag (MHz)', ha='center')
        fig.text(0.05, 0.5, r'Autocorrelation', va='center', rotation='vertical')

        ax2.scatter(10,1200,facecolors='none', edgecolors='none',label=r'\textbf{a}')
        ax2.legend(loc='upper left',handlelength=0, handletextpad=-0.5, bbox_to_anchor=(0.01,0.95 ),frameon=False,markerscale=0, fontsize=8)

        ax4.scatter(10,0.1,facecolors='none', edgecolors='none',label=r'\textbf{b}')
        ax4.legend(loc='upper left',handlelength=0, handletextpad=-0.5, bbox_to_anchor=(0.022,0.95 ),frameon=False,markerscale=0, fontsize=8)

        plt.show()

    print("done")
