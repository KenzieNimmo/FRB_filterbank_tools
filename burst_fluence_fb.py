"""
Filterbank of burst -> fluence and peak flux density of the burst
Kenzie Nimmo 2020
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import optparse
import re 

def radiometer(tsamp, bw, npol, SEFD):
    """
    radiometer(tsamp, bw, npol, Tsys, G):
    tsamp is the time resolution in milliseconds
    bw is the bandwidth in MHz
    npol is the number of polarizations
    Tsys is the system temperature in K (typical value for Effelsberg = 20K)
    G is the telescope gain in K/Jy (typical value for Effelsberg = 1.54K/Jy)
    """

    return (SEFD)*(1/np.sqrt((bw*1.e6)*npol*tsamp*1e-3))


def fluence_flux(arr, bw, t_cent, width, tsamp, SEFD, offpulse):
    """
    fluence_flux(arr, bw, t_cent, width, tsamp, offpulse)
    arr is the burst dynamic spectrum
    bw is the bandwidth in MHz
    t_cent is the peak time of the burst found using the 2D Gaussian fit (or by e\
ye)
    width is the FWHM duration of the burst found using the 2D Gaussian fit (or b\
y eye)
    tsamp is the sampling time of the data in seconds
    offpulse is the pickle file containing the offpulse times

    Idea is to subtract mean and divide by the rms to normalize the time series
    (making the noise ~1 and so the height of the signal is equal to the S/N)
    Then to convert to physical units (Jy ms), we use the radiometer equation.

    Also use same method to determine peak flux in physical units (Jy)
    """

    with open(str(offpulse),'rb') as f:
        offtimes=pickle.load(f)

    f.close()
    #offtimes=offtimes[offtimes<2400]
    tsamp*=1e3 #milliseconds

    conv = 2.355
    width=int((width*2./conv))
    t_cent = int(t_cent)

    profile = np.mean(arr,axis=0)
    spec = np.mean(arr[:,(t_cent-width):(t_cent+width)],axis=1)
    offprof = np.mean(arr[:,offtimes],axis=0)
    offspec = np.mean(arr[:,offtimes],axis=1)
    mean = np.mean(offprof)
    meanspec=np.mean(offspec)
    offprof-=mean
    profile-=mean
    spec-=meanspec
    std = np.std(offprof)
    stdspec=np.std(offspec)
    offprof /=std
    profile/=std
    spec/=stdspec

    profile_burst = profile[(t_cent-width):(t_cent+width)]
    spec_burst = spec
    

    fluence= np.sum(profile_burst*radiometer(tsamp,bw,2,SEFD)*tsamp) #fluence
    flux=np.max(profile_burst*radiometer(tsamp,bw,2, SEFD)) #peak flux density
    prof_flux=profile*radiometer(tsamp,bw,2, SEFD)
    spec_flux=spec_burst*radiometer(tsamp,bw,2, SEFD)
    return fluence, flux, prof_flux, spec_flux

def energy_iso(fluence,distance_lum):
    """                                                                              Following Law et al. (2017)                                                  
    fluence in Jy ms                               
    distance_lum in Mpc                                                          
    At the moment bw not needed, units are therefore erg Hz^-1
    """
    #convert Jy ms to J s                                                        
    fluence_Jys = fluence*1e-3                                                   
    #convert Mpc to cm                                                           
    distance_lum_cm = 3.086e24*distance_lum
    return fluence_Jys*4*np.pi*(distance_lum_cm**2)*1e-23

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [options] infile', \
                description="Fluence and peak flux density of FRB. Input the pickle file o\
utput from fit_burst_fb.py.")
    parser.add_option('-n', '--burst_no', dest='burst_no', type='int', \
                      help="Burst number used as an index.", default=None)
    parser.add_option('-S', '--SEFD', dest='SEFD', type='float', \
                      help="System Equivalent Flux Density [Jy].", default=None)
    parser.add_option('-d', '--distance', dest='distance', type='float', \
                      help="Distance to the FRB for energy calculation in Mpc (not required).", default=None)
    parser.add_option('-u', '--uncorr', dest='uncorr', action="store_true", \
                      help="If -u option is used, use the uncorrected array (othe\
rwise use the masked+bandpass corrected array).", default=False)

    (options, args) = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        sys.exit(1)
    elif len(args)!=1:
        sys.stderr.write("Only one input file must be provided!\n")
    else:
        options.infile = args[-1]

    if options.burst_no is None:
        sys.stderr.write("A burst index must be provided." \
                            "(Use -n/--burst_no on command line).\n")
        sys.exit(1)

    if options.SEFD is None:
        sys.stderr.write("An SEFD must be provided." \
                            "(Use -S/--SEFD on command line).\n")
        sys.exit(1)


    pklfilename = str(options.infile)
    with open(pklfilename, 'rb') as f:
        bursts=pickle.load(f)

    f.close()

    path='./'
    IDs_ordered=[str(options.burst_no)]
    for burst in IDs_ordered:
        burst=str(burst)

        arr_corr = bursts[burst]['array_corrected']
        arr_uncorr = bursts[burst]['array_uncorrected'] #no rfi masking or bp correction
        if options.uncorr == True:
            arr=arr_uncorr
        else:
            arr=arr_corr

        tsamp = bursts[burst]['t_samp']
        freqs = bursts[burst]['freqs']
        nchan=len(freqs)
        fres=np.abs((freqs[-1]-freqs[0])/(nchan-1))
        bw=np.abs(freqs[-1]-freqs[0])+fres
        
        t_cent = bursts[burst]['centre_bin']
        t_fwhm = bursts[burst]['width_bin']
        
        basename=re.search('(.*).pkl',pklfilename).group(1)
        offpulse = "%s_offpulse_time.pkl"%basename

        fluence, flux, prof_flux, spec_flux = fluence_flux(arr, bw, t_cent, t_fwhm, tsamp,options.SEFD, offpulse)


        print("Fluence:", fluence, "Jy ms")
        print("Flux Density:", flux, "Jy")
        if options.distance!=None:
            specenerg = energy_iso(fluence,options.distance)
            print("Spectral energy density:", specenerg, "erg Hz^{-1}")
            
        fits={}
        burst_properties={}
        fits['array_corrected']=bursts[burst]['array_corrected']
        fits['array_uncorrected']=bursts[burst]['array_uncorrected']
        fits['mask']=bursts[burst]['mask']
        fits['t_samp']=bursts[burst]['t_samp']
        fits['freqs']=bursts[burst]['freqs']
        fits['centre_bin']=bursts[burst]['centre_bin']
        fits['width_bin']=np.abs(bursts[burst]['width_bin'])
        if bursts[burst].get('scint_bw')!=None:
            fits['scint_bw']=bursts[burst]['scint_bw']
            fits['std_ACF_fit']=bursts[burst]['std_ACF_fit']
            fits['chisq_ACF_fit']=bursts[burst]['chisq_ACF_fit']
            fits['dof_ACF_fit']=bursts[burst]['dof_ACF_fit']
            fits['pval_ACF_fit']=bursts[burst]['pval_ACF_fit']
            fits['pcov_ACF_fit']=bursts[burst]['pcov_ACF_fit']
            fits['ACF']=bursts[burst]['ACF']
            fits['lorentzian']=bursts[burst]['lorentzian']
            fits['freq_lag']=bursts[burst]['freq_lag']
            fits['freq_lorentz']=bursts[burst]['freq_lorentz']

        fits['fluence']=fluence
        fits['peakfluxdens']=flux
        fits['prof_flux']=prof_flux
        fits['spec_flux']=spec_flux
        if options.distance!=None:
            fits['specenergdens']=specenerg


        burst_properties[burst] = fits
        
            

        with open(pklfilename, 'wb') as f:
            pickle.dump(burst_properties, f)
