"""
Interactive RFI masker using the filterbank file, and outputs the masked array that can be read into other python scripts
Kenzie Nimmo 2020
"""

import numpy as np
import matplotlib.pyplot as plt
import filterbank_to_arr
from pylab import *
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import filterbank
import pickle

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

rows=2
cols=2
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, wspace=0., hspace=0., height_ratios=[0.5,]*(rows-1)+[2,], width_ratios=[5,]+[1,]*(cols-1))
ax1 = plt.subplot(gs[2]) #dynamic spectrum
ax2 = plt.subplot(gs[0],sharex=ax1) #profile
ax3 = plt.subplot(gs[-1],sharey=ax1) #spectrum

fil=filterbank.filterbank('test.fil')
total_N=1000
arr = filterbank_to_arr.filterbank_to_np('test.fil',None,total_N)
spectrum=np.mean(arr,axis=1)
profile=np.mean(arr,axis=0)
freqbins=np.arange(0,arr.shape[0],1)
freqs=fil.frequencies
threshold=np.amax(arr)-(np.abs(np.amax(arr)-np.amin(arr))*0.99) #create a threshold for imaging the brightest points (top 1%) as a different colour.
#bright_points = np.ma.masked_greater(arr, threshold) #wherever this is true are the bright points

cmap = mpl.cm.binary

ax1plot = ax1.imshow(arr,aspect='auto',vmin=np.amin(arr),vmax=threshold,cmap=cmap,origin='lower',interpolation='nearest',picker=True)
cmap.set_over(color='pink')
ax1.set_xlim(0,total_N)
ax2plot, = ax2.plot(profile, 'k-',alpha=1.0,zorder=1)
ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
ax2.tick_params(axis='x', labelbottom='off', top='off')
y_range = profile.max() - profile.min()
ax2.set_ylim(profile.min()-y_range*0.1, profile.max()*1.1)
ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False, bottom=True, top=True, left=True, right=True)
ax3plot, = ax3.plot(spectrum, freqbins, 'k-',zorder=2)
ax3.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='off')
ax3.tick_params(axis='y', labelleft='off')
ax3.set_ylim(freqbins[0], freqbins[-1])
x_range = spectrum.max() - spectrum.min()
ax3.set_xlim(-x_range/4., x_range*6./5.)

fig.add_subplot(ax1)
fig.add_subplot(ax2)
fig.add_subplot(ax3)

begin_chan = []
mask_chan = []
def onclick(event):
    tb = get_current_fig_manager().toolbar
    if tb.mode == '':
        y1= event.ydata
        arr=ax1plot.get_array()
        vmin = np.amin(arr)
        index=find_nearest(freqbins,y1)
        begin_chan.append(index)

def onrel(event,ithres):
    tb = get_current_fig_manager().toolbar
    if tb.mode == '':
        y2= event.ydata
        arr=ax1plot.get_array()
        vmin = np.amin(arr)
        index2=find_nearest(freqbins,y2)
        if begin_chan[-1] > index2:
            arr[index2:begin_chan[-1],:]=vmin-100
        else:
            arr[begin_chan[-1]:index2,:]=vmin-100
        mask = arr<vmin-50
        arr = np.ma.masked_where(mask==True,arr)
        ax1plot.set_data(arr)
        profile = np.mean(arr,axis=0)
        ax2plot.set_data(np.arange(0,total_N,1),profile)
        ax3plot.set_data(np.mean(arr,axis=1),freqbins)
        threshold=np.amax(arr)-(np.abs(np.amax(arr)-np.amin(arr))*ithres)
        ithres-=0.05
        ax1plot.set_clim(vmin=np.amin(arr),vmax=threshold)
        spectrum =  np.mean(arr,axis=1)
        ax3.set_xlim(np.amin(spectrum),np.amax(spectrum))
        y_range = profile.max() - profile.min()
        ax2.set_ylim(profile.min()-y_range*0.1, profile.max()*1.1)

        cmap.set_over(color='pink')
        plt.draw()
        if begin_chan[-1] > index2:
            mask_chan.append([i for i in np.arange(index2,begin_chan[-1],1)])
        else:
            mask_chan.append([i for i in np.arange(begin_chan[-1],index2,1)])


ithres=0.9
plt.connect('button_press_event', onclick)
plt.connect('button_release_event', lambda event: onrel(event, ithres))
plt.show()

mask_chans = [item for sublist in mask_chan for item in sublist]
mask_chans = list(dict.fromkeys(mask_chans))
with open('mask.pkl', 'wb') as f:
    pickle.dump(mask_chans, f)
