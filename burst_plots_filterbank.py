
"""
Kenzie Nimmo 2020
"""
import os
import argparse
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import psrchive
import scipy.ndimage
from matplotlib.ticker import MultipleLocator
import astropy.constants as cc
import scipy.ndimage.filters as filters
import subprocess
import pandas as pd
import psr_utils
import matplotlib.colors
import itertools
import pickle
import sys
import matplotlib.patheffects as PathEffects


mpl.rcParams['font.size'] = 7
#mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['legend.fontsize'] = 7
mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['text.latex.preamble'] = r'\usepackage[T1]{fontenc}'


def plot(pkl,bind,favg,tavg,width,fig,plot_grid_idx,index,threshold,plot_spectrum=False,addmask=None):

    """
    Creates the montage of bursts plot. ds, extent and w are outputs from load_archive. plot_grid_idx gives the plot grid structure with index number of plots/bursts. fig is the plt.figure(). t_cent is the central time of the burst with time fwhm, t_fwhm. f_cent and f_fwhm are the equivalent in frequency. time values are given in number of bins. frequency in MHz.
    plot_spectrum = True will plot the spectrum (amplitude vs frequency).
    smoothed=True will return a smoothed pulse profile.
    If broadband burst: select broadband=True and do not give f_cent or f_fwhm. If narrow band, select broadband=False and give f_cent and f_fwhm from the 2D gaussian fit code.
    """
    with open(pkl,'rb') as f:
        pklfile=pickle.load(f)
    f.close()

    arr_corr=pklfile['array_corrected']
    arr_uncorr=pklfile['array_uncorrected']
    mask = pklfile['mask']
    res_t = pklfile['t_samp']
    freqs = pklfile['freqs']
    nchan_tot = len(freqs)
    mask=nchan_tot-mask
    res_f = (freqs[-1]-freqs[0])/nchan_tot
    fmin = freqs[0]-(res_f/2.)
    fmax = freqs[-1]-(res_f/2.)
    t_cent = pklfile['centre_bin']
    t_fwhm = pklfile['width_bin']
    prof_flux = pklfile['prof_flux']
    spec_flux=pklfile['spec_flux']
    spec_flux[mask]=0
    nbins = np.shape(arr_corr)[1]

    conv = 2.355 #convert FWHM to sigma

    if favg == None:
        favg = int(1)
    if tavg == None:
        tavg = int(1)

    arr=arr_uncorr
    spectrum=np.mean(arr,axis=1)
    if addmask!=None:
        for i in range(len(addmask)):
            mask=np.append(mask,addmask[i])
    #arr[mask]=0
    maskarray=np.ones_like(spectrum)
    maskarray[mask]=0

    #downsample in time
    if tavg>1:
        tavg=float(tavg)
        binstoremove=(nbins/tavg)-np.int(nbins/tavg)
        binstoremove*=tavg
        endind=int(nbins-binstoremove)
        arr_trim=arr[:,:endind]
        prof_flux_trim=prof_flux[:endind]
        newnbins=np.shape(arr_trim)[1]
        tavgbins=newnbins/tavg
        arr=np.array(np.column_stack([np.mean(subint, axis=1) for subint in np.hsplit(arr_trim,tavgbins)]))
        tavg=int(tavg)
        prof_flux =np.mean(prof_flux_trim.reshape(-1, tavg), axis=1) #determining the mean flux for plot ax1

    #downsample in frequency
    if favg>1:
        favg=float(favg)
        if (nchan_tot/favg)-int(nchan_tot/favg)!=0:
            print("The total number of channels is %s, please choose an fscrunch value that divides the total number of channels.")
            sys.exit()
        else:
            newnchan=nchan_tot/favg
            arr=np.array(np.row_stack([np.mean(subint, axis=0) for subint in np.vsplit(arr,newnchan)]))
            favg=int(favg)
            maskarray=np.mean(maskarray.reshape(-1,favg),axis=1)
            spec_flux=np.mean(spec_flux.reshape(-1,favg),axis=1)

    mask=np.where(maskarray==0)[0]
    arr[mask]=0
    #determine the time and frequency resolution of the plot based on downsampling
    res_t*=tavg
    res_f*=favg

    f_l_bin=0
    f_h_bin=nchan_tot-1

    extent=np.zeros(4) #time_min,time_max,freq_min,freq_max
    peak = np.int(np.ceil(t_cent))/tavg #time bin containing peak of burst
    extent[0] = - width / 2. # width is time window around the burst
    extent[1] = width / 2.
    extent[2] = fmin
    extent[3] = fmax
    centre = np.ceil(width / 2. / (res_t*1e3)) #bin number of the centre of window


    #cut out the relevant part of the array (within extent)
    t_fwhm/=tavg
    #spectrum = np.mean(arr[:,int(np.ceil(peak-t_fwhm/2.)):int(np.ceil(peak+t_fwhm/2.))], axis=1)
    spectrum=spec_flux
    #maskind=np.where(spectrum==0)
    #maskspec=spectrum==0
    maskspec=maskarray
    spectrum=np.ma.masked_where(maskspec==False,spectrum)
    arr=arr[:,int(np.ceil(peak-centre)):int(np.ceil(peak+centre))]
    #ts = np.mean(arr,axis=0) #time series
    ts=prof_flux[int(np.ceil(peak-centre)):int(np.ceil(peak+centre))]
    """
    #This part is to centre the burst. Different cases depending on where the burst is, relative to the centre, to begin with.
    t1 = int(peak + centre)
    if (peak-centre) < 0:
        t0 = np.abs(int(peak-centre))
        wrap = ts[(t_bins-t0):]
        full_wrap = full_ts[(t_bins-t0):]
        wrap_ds = ds[fmin_bin : fmax_bin, (t_bins-t0):]
        ds_end = ds[fmin_bin : fmax_bin, 0:t1]
        ts_end = ts[0:t1]
        full_ts_end = full_ts[0:t1]
        ts = np.append(wrap,ts_end)
        full_ts = np.append(full_wrap,full_ts_end)
        ds = np.append(wrap_ds,ds_end,axis=1)

    elif t1 > t_bins:
        t0 = int(peak-centre)
        wrap_num = t1 - ts.shape[0]
        wrap = ts[0:wrap_num]
        full_wrap = full_ts[0:wrap_num]
        wrap_ds = ds[fmin_bin : fmax_bin, 0 : wrap_num]
        ds_beg = ds[fmin_bin : fmax_bin, t0 : ]
        ts_beg = ts[t0:]
        full_ts_beg = full_ts[t0:]
        ts = np.append(ts_beg,wrap)
        full_ts = np.append(full_ts_beg,full_wrap)
        ds = np.append(ds_beg,wrap_ds,axis=1)
    else:
        t0 = int(peak-centre)
        ts = ts[t0:t1]
        full_ts = full_ts[t0:t1]
        ds = ds[fmin_bin : fmax_bin, t0 : t1]
    """
    # creating the figure structure
    rows=2
    cols=1
    if plot_spectrum: cols += 1
    plot_grid = gridspec.GridSpecFromSubplotSpec(rows, cols, plot_grid_idx, wspace=0., hspace=0.,\
                                                 height_ratios=[1,]*(rows-1)+[2,], width_ratios=[5,]+[1,]*(cols-1))

    ax1 = plt.Subplot(fig, plot_grid[rows-1,0])
    ax2 = plt.Subplot(fig, plot_grid[rows-2,0], sharex=ax1)
    if plot_spectrum: ax3 = plt.Subplot(fig, plot_grid[rows-1,1], sharey=ax1)
    units = ("GHz", "ms")

    # plotting the waterfall plot with the burst at the centre, arr
    cm1 = mpl.colors.ListedColormap(['black','red'])
    cm2 = mpl.colors.ListedColormap(['black','blue'])

    vmin = np.amin(arr) #min of the array
    vmax=np.amax(arr)   #max of the array
    zapped = mask #identifying the frequencies that were masked
    cmap = plt.cm.gist_yarg
    cmap.set_bad((1, 0, 0, 1)) #set color for masked values
    zap_size = int(arr.shape[1]/18)
    arr[zapped,:zap_size] = vmin-1000
    mask1=arr<vmin-600
    mask1 = np.ma.masked_where(mask1==False,mask1)
    cmap.set_over(color='white') #set color used for high out-of-range values (the zapped values have been set to NaN)
    arr[zapped,zap_size:] = vmax+1000.

    ax1.imshow(arr, cmap=cmap, origin='lower', aspect='auto', interpolation='nearest', vmin=vmin, vmax=vmax, extent=extent)
    ax1.imshow(mask1, cmap=cm1, origin='lower', aspect='auto', interpolation='nearest', vmin=0, vmax=1, extent=extent)

    if width: ax1.set_xlim(-width/2.-0.001, width/2.+0.001)
    ax1.set_ylim(fmin, fmax)


    #Label only edge plots
    if index % ncols == 0:
        ax1.set_ylabel(r'${\rm Frequency}\ ({\rm GHz})$')#.format(units[0]))
        yticks=np.array([2.2,2.22,2.24,2.26,2.28,2.30])
        ax1.set_yticks(yticks*1000.)
        ax1.set_yticklabels([r'$2.20$',r'$2.22$',r'$2.24$',r'$2.26$',r'$2.28$',r'$2.30$'])
        ax2.set_ylabel(r'${\rm Flux\ Density}\ ({\rm Jy})$')
    else:ax1.tick_params(axis='y', labelleft='off')


    if (index <threshold) and width: ax1.tick_params(axis='x', labelbottom='off')
    else:
        ax1.set_xlabel(r'${\rm Time}\ ({\rm ms})$')#.format(units[1]))
        xticks=([(-width/2.),(-width/4.),0,(width/4.),(width/2.)])
        ax1.set_xticks(xticks)
        #ax1.set_xticklabels()
    ax1.xaxis.set_minor_locator(MultipleLocator(0.5))

    #plot pulse profile (time series)
    x = np.linspace(extent[0], extent[1], ts.shape[0])
    ax2.plot(x, ts, 'k-',alpha=1.0,zorder=1)
    ax2.set_yticks(np.array([0,round(ts.max(),2)]))
    box = ax2.get_position()
    #ax2.set_position([box.x0, box.y0, box.width*0.65, box.height])
    point = ax2.scatter(x[0], ts[0], facecolors='none', edgecolors='none')
    legend_x = 0.82
    legend_y = 0.72
    #ax2.scatter(10,1200,facecolors='none', edgecolors='none',label=(r'${\rm B%s}$')%bind)#.format(burst_n)))
    #ax2.scatter(10,1200,facecolors='none', edgecolors='none',label=(r'$%.3f\ {\rm MHz}$')%res_f)#(r'${0:.3f}\ {\rm MHz}$'.format(f)))
    #ax2.scatter(10,1200,facecolors='none', edgecolors='none',label=(r'$%.3f\ {\rm ms}$')%(res_t*1e3))#(r'${0:.3f}\ {\rm ms}$'.format(t)))
    txtfres=ax2.text(0.9,0.85,r'$%.3f\ {\rm MHz}$'%res_f, ha='center', va='center', transform = ax2.transAxes)
    txtfres.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
    txttres=ax2.text(0.9,0.68,r'$%.3f\ {\rm ms}$'%(res_t*1e3), ha='center', va='center', transform = ax2.transAxes)
    txttres.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
    #ax2.tick_params(axis='y', which='both', right='off')
    ax2.tick_params(axis='x', labelbottom='off', top='off')
    y_range = ts.max() - ts.min()
    ax2.set_ylim(-y_range/3., ts.max()*1.1)
    #ax2.hlines(y=(-y_range/3.),xmin=-1,xmax=1,lw=10,color='#d6ffff',zorder=0.8)
    ax2.hlines(y=-y_range/3.,xmin=(-t_fwhm*2./conv*res_t*1e3),xmax=(t_fwhm*2./conv*res_t*1e3), lw=10,color='#48D1CC',zorder=0.8,alpha=0.3)#, alpha=0.3)#color='#FFD700',zorder=3)
    ax2.hlines(y=-y_range/3., xmin=(-t_fwhm/2.*res_t*1e3),xmax=(t_fwhm/2.*res_t*1e3), lw=10, color='#48D1CC',zorder=0.9)

    legend_x = 0.76
    legend_y = 0.74
    #legend2=ax2.legend(loc='center left',handlelength=0, handletextpad=-0.5, bbox_to_anchor=(legend_x, legend_y),frameon=False,markerscale=0)
    legend1=ax2.legend((point,point), ((r'${\rm B%s}$')%bind,""), loc='upper left',handlelength=0,bbox_to_anchor=(0.02, 0.95), handletextpad=-0.5,frameon=False,markerscale=0)
    #ax2.add_artist(legend2)
    ax2.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False, bottom=True, top=True, left=True, right=True)

    #plot spectrum (amplitude vs freq) only if plot_spectrum=True
    if plot_spectrum:
        y = np.linspace(extent[2], extent[3], spectrum.size)
        ax3.plot(spectrum, y, 'k-',zorder=2)
        ax3.tick_params(axis='x', which='both', top='off', bottom='on', labelbottom='on')
        ax3.tick_params(axis='y', labelleft='off')
        ax3.set_ylim(fmin, fmax)
        x_range = spectrum.max() - spectrum.min()
        ax3.set_xticks(np.array([0,round(np.max(spectrum),2)]))
        ax3.set_xticklabels([r'$0.00$',r'$%.2f$'%np.max(spectrum)])
        if (index <threshold) and width:
            ax3.set_xlim(-x_range/3., x_range*6./5.)
        else:
            ax3.tick_params(axis='x', pad=6.0, length=6.0, direction='out')
            ax3.set_xlabel(r'${\rm Flux\ Density}\ ({\rm Jy})$')
            ax3.set_xlim(-x_range/2., x_range*6./5.)

    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    if plot_spectrum: fig.add_subplot(ax3)
    return



if __name__ == '__main__':
    #adapt the figure parameters as needed:
    nrows = 3
    ncols = 2 # change this to suit the number of bursts you want to plot
    threshold=(nrows*ncols)-ncols #for axes on plots
    width = 20. #Time window around the burst in ms.
    plot_grid = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1) #grid of burst plots
    fig = plt.figure(figsize=[7,7]) #defines the size of the plot
    archive_path = '/data2/nimmo/FRB_R3_2019/20190619/FRBR3/burst_archives/final_psr_archives/' #directory containing burst archives (create using dspsr)
    IDs_ordered = [1,2,3,4,5,6]
    bursts = {"1": {'pkl': 'burst_411_sband_nomask.pkl', 'favg':2, 'tavg':None, 'addmask':[117,58,59,106,105,104,103]},
              "2": {'pkl': 'burst_12118_sband_nomask.pkl', 'favg':2, 'tavg':2, 'addmask':[117]},
              "3": {'pkl': 'burst_607_sband_nomask.pkl', 'favg':2, 'tavg':None, 'addmask':[117]},
              "4": {'pkl': 'burst_5888_sband_nomask.pkl', 'favg':2, 'tavg':2, 'addmask':[117]},
              "5": {'pkl': 'burst_623_sband_nomask.pkl', 'favg':2, 'tavg':2, 'addmask':[117]},
              "6": {'pkl': 'burst_872_sband_nomask.pkl', 'favg':2, 'tavg':None, 'addmask':[117]}

}


    idx = 0
    for burst in IDs_ordered:
        burst=str(burst)
        plot(bursts[burst]['pkl'],burst,bursts[burst]['favg'],bursts[burst]['tavg'],width,fig,plot_grid[idx],idx,threshold,plot_spectrum=True,addmask=bursts[burst]['addmask'])


        idx+=1


#plot_dispersed(archive_name,t_cent,deltat,index,fig,plot_grid_idx,y_axis=False)
    fig.subplots_adjust(hspace=0.08, wspace=0.05, left=0.09,right=.96,bottom=.1,top=.97)
    #plt.show()
    fig.savefig("bursts_KN.pdf", format = 'pdf', dpi=300)
