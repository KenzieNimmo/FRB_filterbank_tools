# FRB_filterbank_tools
A selection of python scripts used to manipulate filterbank files in the analysis of FRB data

extract.py -o output_filename -b begin_sample_number -d duration_number_of_samples input_filename: extract a filterbank file from input_filename called output_filename, beginning at sample number begin_sample_number for a total number of samples duration_number_of_samples.
<<May be an error with the start time in the new filterbank header information -- look into this>>

RFI_offpulse.py: is an interactive tool used to manually mask RFI in your data, and define where the burst is (so you know an offpulse region for bandpass correction).
This outputs a pickle file (*basename*.pkl) that contains an array of the data unmasked and without bandpass correction, a copy of the mask, an array of the data masked and bandpass corrected, the sample time of your data and a list of frequencies. 
To see input options: >> python RFI_offpulse.py -h 
To use: click r to refresh if dynamic range shows mostly pink on plot
        click x and click and drag your mouse to define where the burst is in the time series
        click y to remove this if you make a mistake or change your region
        to zap RFI, click on the dynamic spectrum and drag to include all of the channels you want to remove. Equivalently you can click and drag on the spectrum to remove also.


python fit_bursts_fb.py *basename*.pkl : fits a 2D Gaussian fit (freq and time) to your burst. There is an option to provide a guess for the initial fit. There are diagnostic plots for you to convince yourself that the fit is good. This script uses fit_2D_Gaussian.py, fit_repeater_bursts.py and plotfuncs.py. Output is a pickle file (*basename*.pkl) containing an array of the data unmasked and without bandpass correction, a copy of the mask, an array of the data masked and bandpass corrected, the sample time of your data, a list of frequencies, the centre time bin of the fit, and the FWHM duration of the fit.
To see input options: >>python fit_bursts_fb.py -h
<<I have implemented barycentre correcting the peak time but it requires python 3 and astropy version > 2.2 which is not yet working on dragnet, therefore this is currently untested.>>


python ACF_filterbank.py *basename*.pkl : autocorrelation function analysis of the spectrum of your burst, using the peak time and width determined using fit_bursts_fb.py. Output is a pickle file (*basename*.pkl) containing an array of the data unmasked and without bandpass correction, a copy of the mask, an array of the data masked and bandpass corrected, the sample time of your data, a list of frequencies, the centre time bin of the fit, the FWHM duration of the fit, the scintillation bandwidth, the standard deviation of the ACF, the chisq of the Lorentzian fit, the degrees of freedom of the Lorentzian fit, the pvalue of the Lorentzian fit, and the corrected pcov output from curvefit (used for uncertainties).
To see input options: >>python ACF_filterbank.py -h
<<Error propagation in progress>>

python burst_fluence_flux.py *basename*.pkl : calculates the burst fluence and peak flux density using the radiometer equation. Outputs the fluence, peak flux density, profile of the burst in units of Jy, and energy (if you give a distance as input).
To see input options: >>python burst_fluence_flux.py -h

Note you can run burst_fluence_flux.py before or after ACF_filterbank.py in the analysis (they are not dependent on one another). However, before either burst_fluence_flux.py or ACF_filterbank.py you must run RFI_offpulse.py followed by fit_bursts_fb.py.


Final output pickle file *basename*.pkl will contain a dictionary with keys:
['array_corrected'],['array_uncorrected'],['mask'],['t_samp'],['freqs'],['centre_bin'],['width_bin'],['scint_bw'],['std_ACF_fit'],['chisq_ACF_fit'],['dof_ACF_fit'],['pval_ACF_fit'],['pcov_ACF_fit'],['fluence'],['peakfluxdens'],['prof_flux'],['specenergdens']

