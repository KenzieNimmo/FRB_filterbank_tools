"""
python ACF_filterbank.py *basename*.pkl 
KN 2020
"""
import sys
sys.path.insert(1,'~/FRB_filterbank_tools')
import os
import filterbank
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import rfifind
from scipy.optimize import curve_fit, leastsq
import filterbank_to_arr
import pickle
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

def error(x,rms,nchan,v=None):
    if v is None:
        v= np.ones_like(x)
    x = x.copy()
    x[v!=0] -= x[v!=0].mean()
    errors = np.zeros_like(x)
    for i in range(len(x)):
        if i>0:
            m=shift(v,0,nchan)*shift(v,i,nchan)
            sigmaxsq=(shift(x, 0, nchan)**2)*m*(shift(rms,i,nchan)**2)+(shift(x, i,nchan)**2)*m*(shift(rms,0,nchan)**2)
            sigmaysq=np.sum(sigmaxsq) #this is the square of the error on the top line of the ACF formula
            print(sigmaysq)
            sigmaasq=np.sum(2*(shift(x, 0, nchan)**2)*m*(shift(rms,0,nchan)**2))
            sigmabsq=np.sum(2*(shift(x, i, nchan)**2)*m*(shift(rms,i,nchan)**2))
            sigmar = (0.5*np.sqrt(sigmaasq))/np.sqrt(np.sum((shift(x, 0, nchan)**2)))
            sigmas = (0.5*np.sqrt(sigmabsq))/np.sqrt(np.sum((shift(x, i, nchan)**2)))
            print(sigmar)
            print(sigmas)
            #determine the final error
            y=np.sum(shift(x,0,nchan)*shift(x, i,nchan)*m)
            r=np.sqrt(np.sum(shift(x, 0, nchan)**2*m))
            s=np.sqrt(np.sum(shift(x, i, nchan)**2*m))
            finalACF=y/(r*s)
            sigmafinal= finalACF*np.sqrt((sigmaysq/y**2)+(sigmar**2/r**2)+(sigmas**2/s**2))
            
            print(np.abs(sigmafinal))
            exit()
            errors[i-1]=np.abs(sigmafinal)
    return errors
                                                     



def scint_bw(arr,mask,freqs,burst_peak,burst_width,normalise=False,chan_min=None, chan_max=None):

    nchan_tot = arr.shape[0]
    begin_burst = np.int(burst_peak-(0.5*burst_width))
    end_burst = np.int(burst_peak+(0.5*burst_width))
    if normalise!=False:
            f=open(normalise,'rb')
            offpulse_times=pickle.load(f)
            f.close()
            offspec = np.mean(arr[:,offpulse_times],axis=1)
            offspec_rms = np.zeros_like(offspec)
            for i in range(nchan_tot):
                    #offspec_rms[i]=np.sqrt(((arr[i,offpulse_times]-np.mean(arr[i,offpulse_times]))**2).mean())
                    #offspec_rms[i]=np.std(arr[i,offpulse_times]-np.mean(arr[i,offpulse_times]))

                    #let's try std of on pulse
                    offspec_rms[i]=np.std(arr[i,begin_burst:end_burst]-np.mean(arr[i,begin_burst:end_burst]))


            if chan_min!=None and chan_max!=None:
                    offspec = np.mean(arr[chan_min:chan_max,offpulse_times],axis=1)
 

    spec = np.mean(arr[:,begin_burst:end_burst],axis=1)
    if chan_min!=None and chan_max!=None:
            spec=np.mean(arr[chan_min:chan_max,begin_burst:end_burst],axis=1)

    #spec = spec.data

    if np.all(mask)!=None:
            spec[mask]=0
    
    if chan_min!=None and chan_max!=None:
            nchan=chan_max-chan_min
    else:
            nchan=nchan_tot
    channel_list = np.arange(0,nchan,1)
    f_res = np.abs( (freqs[-1]-freqs[0])/nchan_tot) #MHz
    delta_f=np.linspace(0,nchan,nchan)*f_res
    if chan_min!=None and chan_max!=None:
            delta_f=np.arange(0,chan_max-chan_min,1)*f_res
    #print(delta_f)
    good = np.ones(len(spec))
    if np.all(mask)!=None:
            good[mask]=0 #good is an array of 1 and 0 where 1 is good and 0 is bad (so we ignore the bad)
    ACF=autocorr(spec, nchan=nchan, v=good)
    ACFerror=error(spec,offspec_rms,nchan=nchan,v=good)
    
    rev = ACF[::-1]
    ACF_negtopos = np.append(rev,ACF)
    plt.plot(delta_f,ACF+1.1, 'r', label='ACF', drawstyle='steps-mid')
    if normalise!=False:
            ACFoff = autocorr(offspec, nchan=nchan, v=good)
            plt.plot(delta_f,ACFerror, 'b', label='ACFerror',drawstyle='steps-mid')
    #plt.xlim((0,10))
    plt.ylabel('ACF')
    plt.xlabel('Freq lag (MHz)')

    plt.legend()
    plt.show()

    if normalise!=False: return ACF, delta_f, ACFerror
    else: return ACF, delta_f



def lorentzian_fit(ACF, delta_f, nbins_cent,ACFoff=None):
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

    #define the lorentzian function
    def lorentz(x,gamma,y0, c):
        return (y0*gamma)/(((x)**2)+gamma**2)+c

    #uncert_region = ACF[70:200]
    #std = np.std(uncert_region)
    if np.all(ACFoff)!=None:
            errors = ACFoff[0:nbins_cent]
            reversed_arr = errors[::-1]
            std_array = np.append(errors,reversed_arr)

    #std_array = np.zeros_like(ACFc)
    #std_array+=std
    
    popt, pcov = curve_fit(lorentz, delta_fc, ACFc,sigma=std_array)#, sigma=std_array,maxfev=2000) #note we set the sigma to be the standard deviation calculated above, this is an outdated version of curve_fit but in a newer version we could set absolute_sigma = True to avoid curve_fit scaling the errors
    #we have to deal with this mathematically.
    
    chisq = np.sum((lorentz(delta_fc,*popt)-ACFc)**2/(std_array**2))
    dof = len(ACFc)-3 #predicted_chisq
    pval = 1 - stats.chi2.cdf(chisq, dof)
    pval2 = stats.chi2.sf(chisq,dof)
    correct_pcov = pcov*(dof/chisq)
    print "Scintillation BW:", popt[0], "MHz"
    print "std", std_array
    print "chisq", chisq
    print "dof", dof
    print "pval", pval2
    print "pcov", correct_pcov
    
    return ACFc, delta_fc, lorentz(delta_fc,*popt), delta_fc, popt[0], std_array, chisq, dof, pval, correct_pcov

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
        IDs_ordered=[str(sys.argv[2])]
        for burst in IDs_ordered:
                burst=str(burst)
                                
                arr_corr = bursts[burst]['array_corrected']
                arr_uncorr = bursts[burst]['array_uncorrected'] #no rfi masking or bp correction
                mask = bursts[burst]['mask']
                #mask=np.array([np.arange(65,129,1),np.arange(194,258,1),np.arange(323,387,1),np.arange(452,516,1),np.arange(581,645,1)])
                #mask=np.array([np.arange(257,512,1),np.arange(769,1024,1),np.arange(1281,1536,1),np.arange(1793,2048,1),np.arange(2305,2560,1)])
                mask=np.hstack(mask)
                mask=list(mask)
                tsamp = bursts[burst]['t_samp']
                freqs = bursts[burst]['freqs']
                t_cent = bursts[burst]['centre_bin']
                t_fwhm = bursts[burst]['width_bin']
                normalisefile = "%s_offpulse_time.pkl"%basename

                ACF, delta_f, ACFoff = scint_bw(arr_corr,mask,freqs,t_cent,t_fwhm,normalise=normalisefile)

                with open('offpulse_ACF_b6.pkl','wb') as f:
                        pickle.dump(ACFoff,f)
                
                f.close()
                """
                ACF1, delta_f1, ACFoff1 = scint_bw(arr_corr,None,freqs,t_cent,t_fwhm,normalise='b6_I_31khz_zoom_offpulse_time.pkl',chan_min=0,chan_max=256)
                ACF2, delta_f2, ACFoff2 = scint_bw(arr_corr,None,freqs,t_cent,t_fwhm,normalise='b6_I_31khz_zoom_offpulse_time.pkl',chan_min=512,chan_max=768)
                ACF3, delta_f3, ACFoff3 = scint_bw(arr_corr,None,freqs,t_cent,t_fwhm,normalise='b6_I_31khz_zoom_offpulse_time.pkl',chan_min=1024,chan_max=1280)
                ACF4, delta_f4, ACFoff4 = scint_bw(arr_corr,None,freqs,t_cent,t_fwhm,normalise='b6_I_31khz_zoom_offpulse_time.pkl',chan_min=1536,chan_max=1792)
                ACF5, delta_f5, ACFoff5 = scint_bw(arr_corr,None,freqs,t_cent,t_fwhm,normalise='b6_I_31khz_zoom_offpulse_time.pkl',chan_min=2048,chan_max=2304)
                ACF6, delta_f6, ACFoff6 = scint_bw(arr_corr,None,freqs,t_cent,t_fwhm,normalise='b6_I_31khz_zoom_offpulse_time.pkl',chan_min=2560,chan_max=2816)
                
                fig, (ax1,ax2,ax3,ax4,ax5,ax6)=plt.subplots(6, sharex=True)
                ax1.plot(delta_f1,ACF1+1.1, 'r', label='ACF', drawstyle='steps-mid')
                ax1.plot(delta_f1,ACFoff1, 'b', label='ACFoff',drawstyle='steps-mid')
                ax2.plot(delta_f2,ACF2+1.1, 'r', label='ACF', drawstyle='steps-mid')
                ax2.plot(delta_f2,ACFoff2, 'b', label='ACFoff',drawstyle='steps-mid')
                ax3.plot(delta_f3,ACF3+1.1, 'r', label='ACF', drawstyle='steps-mid')
                ax3.plot(delta_f3,ACFoff3, 'b', label='ACFoff',drawstyle='steps-mid')
                ax4.plot(delta_f4,ACF4+1.1, 'r', label='ACF', drawstyle='steps-mid')
                ax4.plot(delta_f4,ACFoff4, 'b', label='ACFoff',drawstyle='steps-mid')
                ax5.plot(delta_f5,ACF5+1.1, 'r', label='ACF', drawstyle='steps-mid')
                ax5.plot(delta_f5,ACFoff5, 'b', label='ACFoff',drawstyle='steps-mid')
                ax6.plot(delta_f6,ACF6+1.1, 'r', label='ACF', drawstyle='steps-mid')
                ax6.plot(delta_f6,ACFoff6, 'b', label='ACFoff',drawstyle='steps-mid')
                

                f=open('lcp_ACF_subband.pkl','wb')
                pickle.dump([ACF1,ACF2,ACF3,ACF4,ACF5,ACF6,ACFoff1,ACFoff2,ACFoff3,ACFoff4,ACFoff5,ACFoff6],f)
                f.close()
                #fig.ylabel('ACF')
                #fig.xlabel('Freq lag (MHz)')
                #fig.text(0.5,0.95,'lcp (B6)',ha='center', va='center')
                plt.show()
                """
                reversed_arr = delta_f[::-1]
                totdelta_f = np.append(-1*reversed_arr,delta_f[0:])
                ACFc, deltafc, lorentz, deltafnew,bw,std, chisq, dof, pval, correct_pcov=lorentzian_fit(ACF, delta_f,35,ACFoff=ACFoff)
                fig = plt.figure(figsize=[183*mm_to_in,100*mm_to_in])
                grid = plt.GridSpec(2, 1, hspace=0.3, wspace=0.3)

                ax2 = fig.add_subplot(grid[0:1,:])
                ax2.plot(totdelta_f,np.append(ACF[::-1],ACF[0:]), color='darkorange', drawstyle='steps-mid')
                #ax2.set_xticks([-128, -64, 0, 64, 128])
                #ax2.set_yticks([0,0.1])
                #ax2.set_yticklabels([r"$0$", r"$0.1$"])
                ax2.set_ylim([np.nanmin(ACF)*1.1,np.nanmax(ACF)*1.1])
                ax4 = fig.add_subplot(grid[1:,:])
                ax4.plot(deltafc, ACFc, color='darkorange',drawstyle='steps-mid')
                ax4.plot(deltafnew, lorentz, color='darkgreen')
                ax4.axvline(x=bw,color='k',dashes=(5,2), lw=1.5 )
                #ax4.set_xticks([-1,0,1])
                #ax4.set_yticks([0,0.06])
                #ax4.set_yticklabels([r"$0$", r"$0.06$"])
                ax4.tick_params(axis='y', labelleft=True)
                ax4.set_xlim([-1.5,1.5])
                ax4.set_ylim([np.amin(ACFc)*1.1,np.amax(ACFc)*1.1])
                
                fig.text(0.5, 0.03, r'Frequency Lag (MHz)', ha='center')
                fig.text(0.05, 0.5, r'Autocorrelation', va='center', rotation='vertical')
                
                ax2.scatter(10,1200,facecolors='none', edgecolors='none',label=r'\textbf{a}')
                ax2.legend(loc='upper left',handlelength=0, handletextpad=-0.5, bbox_to_anchor=(0.01,0.95 ),frameon=False,markerscale=0, fontsize=8)
                
                ax4.scatter(10,0.1,facecolors='none', edgecolors='none',label=r'\textbf{b}')
                ax4.legend(loc='upper left',handlelength=0, handletextpad=-0.5, bbox_to_anchor=(0.022,0.95 ),frameon=False,markerscale=0, fontsize=8)

                plt.show()
                plt.savefig('ACF_%s.pdf'%basename, format='pdf')

        

                fits={}
                burst_properties={}
                fits['array_corrected']=bursts[burst]['array_corrected']
                fits['array_uncorrected']=bursts[burst]['array_uncorrected']   
                fits['mask']=bursts[burst]['mask']
                fits['t_samp']=bursts[burst]['t_samp']
                fits['freqs']=bursts[burst]['freqs']
                fits['centre_bin']=bursts[burst]['centre_bin']
                fits['width_bin']=np.abs(bursts[burst]['width_bin'])
                fits['scint_bw']=bw
                fits['std_ACF_fit']=std
                fits['chisq_ACF_fit']=chisq
                fits['dof_ACF_fit']=dof
                fits['pval_ACF_fit']=pval
                fits['pcov_ACF_fit']=correct_pcov
                fits['ACF']=ACF
                fits['lorentzian']=lorentz
                fits['freq_lag']=delta_f
                fits['freq_lorentz']=deltafc

                if bursts[burst].get('fluence'):
                        fits['fluence']=bursts[burst]['fluence']
                        fits['peakfluxdens']=bursts[burst]['peakfluxdens']
                        fits['prof_flux']=bursts[burst]['prof_flux']
                        if bursts[burst].get('specenergdens'): fits['specenergdens']=bursts[burst]['specenergdens']

                burst_properties[burst] = fits

        f=open("%s.pkl"%basename, "wb")
        pickle.dump(burst_properties,f)
        f.close()
        print("done")
