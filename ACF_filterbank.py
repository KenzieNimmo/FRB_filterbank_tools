"""
python ACF_filterbank.py *basename*.pkl 
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
import sys
import re 
from scipy import stats
mpl.rcParams['font.size'] = 7
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True

mm_to_in = 0.0393701

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

def scint_bw(arr,mask,freqs,burst_peak,burst_width):

    begin_burst = np.int(burst_peak-(0.5*burst_width))
    end_burst = np.int(burst_peak+(0.5*burst_width))

    spec = np.mean(arr[:,begin_burst:end_burst],axis=1)
    spec = spec.data
    spec[mask]=0

    nchan = arr.shape[0]
    channel_list = np.arange(0,nchan,1)
    f_res = (freqs[0]-freqs[-1])/nchan #MHz
    delta_f=np.linspace(0,nchan,nchan)*f_res

    good = np.ones(len(spec))
    good[mask]=0 #good is an array of 1 and 0 where 1 is good and 0 is bad (so we ignore the bad)

    ACF=autocorr(spec, nchan=nchan, v=good)
    plt.plot(ACF, 'r', label='ACF', drawstyle='steps-mid')

    return ACF, delta_f

def lorentzian_fit(ACF, delta_f, nbins_cent):
    """
    Fits a Lorentzian to the central nbins_cent bins of the ACF
    """
    # only focus on the central part of the ACF
    ACF_cent=ACF[0:nbins_cent]
    delta_cent = delta_f[0:nbins_cent]
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

    popt, pcov = curve_fit(lorentz, delta_fc, ACFc, p0=[0.01,0.2,0], sigma=std_array) #note we set the sigma to be the standard deviation calculated above, this is an outdated version of curve_fit but in a newer version we could set absolute_sigma = True to avoid curve_fit scaling the errors
    #we have to deal with this mathematically.
    
    chisq = np.sum((lorentz(delta_fc,*popt)-ACFc)**2)/(std**2)
    dof = len(ACFc)-3 #predicted_chisq
    pval = 1 - stats.chi2.cdf(chisq, dof)
    correct_pcov = pcov*(dof/chisq)
    print "Scintillation BW:", popt[0], "MHz"

    return ACFc, delta_fc, lorentz(delta_fnew,*popt), delta_fnew, popt[0], std, chisq, dof, pval, correct_pcov

if __name__ == '__main__':
        pklfilename = sys.argv[1] #*.pkl
        basename = re.search('(.*).pkl',pklfilename).group(1)
        #open the dictionary created in fit_bursts_fb.py
        with open(pklfilename, 'rb') as f:
                bursts = pickle.load(f)

        """
        bursts['array_corrected']
        bursts['array_uncorrected']                                                                                                              
        bursts['mask']                                                                                                                           
        bursts['t_samp']                                                                                                                 
        bursts['freqs']                                                                                                                        
        bursts['centre_bin']
        bursts['width_bin']
        """
    
        path='./'
        IDs_ordered=['4']
        for burst in IDs_ordered:
                burst=str(burst)
                                
                arr_corr = bursts[burst]['array_corrected']
                arr_uncorr = bursts[burst]['array_uncorrected'] #no rfi masking or bp correction
                mask = bursts[burst]['mask']
                tsamp = bursts[burst]['t_samp']
                freqs = np.flip(bursts[burst]['freqs'])
                t_cent = bursts[burst]['centre_bin']
                t_fwhm = bursts[burst]['width_bin']

                ACF, delta_f = scint_bw(arr_corr,mask,freqs,t_cent,t_fwhm)
                reversed_arr = delta_f[::-1]
                totdelta_f = np.append(-1*reversed_arr,delta_f[0:])
                ACFc, deltafc, lorentz, deltafnew,bw,std, chisq, dof, pval, correct_pcov=lorentzian_fit(ACF, delta_f,60)

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

                """
                add a final step to add the scintillation bw and all the uncert to *basename*.pkl
                ~std, chisq,\
 dof, pval, correct_pcov~
                """

        

                fits={}
                burst_properties={}
                fits['array_corrected']=bursts[burst]['array_corrected']
                fits['array_uncorrected']=bursts[burst]['array_uncorrected']   
                fits['mask']=bursts[burst]['mask']
                fits['t_samp']=bursts[burst]['t_samp']
                fits['freqs']=bursts[burst]['freqs']
                fits['centre_bin']=bursts[burst]['centre_bin']
                fits['width_bin']=bursts[burst]['width_bin']
                fits['scint_bw']=bw
                fits['std_ACF_fit']=std
                fits['chisq_ACF_fit']=chisq
                fits['dof_ACF_fit']=dof
                fits['pval_ACF_fit']=pval
                fits['pcov_ACF_fit']=correct_pcov

                burst_properties[burst] = fits

        f=open("%s.pkl"%basename, "wb")
        pickle.dump(burst_properties,f)
        f.close()
        print("done")
