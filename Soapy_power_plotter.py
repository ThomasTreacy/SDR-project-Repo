#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:25:33 2022

@author: treacyth

Originally created 4 10 22 to plot two power spectra from
a txt or csv file obtained using soapy_power. 
This is a revision started on 18 10 22.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import argparse
    
def peak_finder(power,freq):
    """
    Uses scipy.signal.find_peaks to find peaks in power array
    freq: array of frequency values corresponding to power array
    """
    peak_indices = find_peaks(power,prominence=1)[0]
    #getting PSD values and frequencies of peaks:
    peak_values = np.array(power[peak_indices])
    peak_freqs = np.array(freq[peak_indices])
    
    #getting 15 highest peaks, so graphs aren't too cluttered
    indices = np.argsort(peak_values)[-15:]
    peak_values = peak_values[indices]
    peak_freqs = peak_freqs[indices]

    return np.array(peak_freqs),np.array(peak_values)

def plot_peaks(peak_freqs,peaks_values):
    """
    plots text labels to peaks at coords in arrays peak_freqs and peak_values
    """
    x = peak_freqs
    y = peaks_values
    for i in range(len(x)):
        #freqs in x are in Hz - converting to Mhz for labelling
        x_label = str(np.round(x[i]/1e6,2)*u.MHz)
        plt.text(x[i],y[i],x_label + "\n" + str(np.round(y[i],2)))
    
def make_data_array(filename, nth_power_spec = 0,file_format = "rtl_power"):
    """
    This function takes a csv file "filename" containing either power spectrum or spectrogram, 
    extracts parameters and puts ONE power spectrum in numpy array
    Also returns upper and lower freq. limits
    "nth_power_spec": plot the nth power spectrum if filename contains spectrogram
    """
    if file_format == "rtl_power":
        spec_dataframe = pd.read_csv(filename,header = None)
        shape = spec_dataframe.shape
        
        #getting parameters and data
        date = spec_dataframe.iloc[0,0]
        timestamp = spec_dataframe.iloc[0,1] 
        f_lower = spec_dataframe.iloc[nth_power_spec,2]
        f_upper = spec_dataframe.iloc[nth_power_spec,3]
        bin_size = spec_dataframe.iloc[0,4]
        samples = spec_dataframe.iloc[0,5]
        data = spec_dataframe.iloc[nth_power_spec,6:shape[1]]
        #data: power spectral density (PSD) values
    
    print("\n")
    print(filename)
    print("Timestamp:", date + timestamp)
    print("freq. range:",f_lower,f_upper)
    print("Bin size:",bin_size/1e3,"kHz")
    print("Samples:",samples)
    print("\n")
    
    #storing data in numpy array
    data = np.array(data)
    return data, f_lower, f_upper

def poly_order_3(x,c,d,e,f):
    return c*x**3 + d*x**2 + e*x + f

def remove_dcoffset(power,freq,func="cubic"):
    """
    This function returns array that is the input array with the central values replaced,
    to remove DC spike.
    power: array of PSD values
    freq: array of freq values corresponding to "power"
    """
    num = 20#size of central interval to remove
    #removing DC spike
    #defining interval :
    lower = int(len(power)/2 - num)
    upper = int(len(power)/2 + num)
    
    #interpolating across gap:
    num = 200#num of points either side of gap to perform fitting on
    #getting data points either side of DC spike gap
    power2 = power[lower - num:lower]
    power3 = power[upper:upper + num]
    freq2 = freq[lower - num:lower]
    freq3 = freq[upper:upper + num]
    power4 = np.concatenate([power2,power3])
    freq4 = np.concatenate([freq2,freq3])
    if func == "linear":
        power[lower:upper] = np.average(power4)
    else:
        pars = curve_fit(poly_order_3,freq4,power4)[0]
        gap = poly_order_3(freq[lower:upper],pars[0],pars[1],pars[2],pars[3])  
        power[lower:upper] = gap
    
    return power, freq

def truncate(arr):
    """
    removes values from start and end of input array, 
    for removing side artefacts from power spectra
    """
    #removing side artefacts
    num = 10
    lower = num
    upper = int(len(arr) - num)
    arr = arr[lower:upper] 

    return arr

def plot_power_spec(paths,peakfinder = True,norm = False,nth_power_spec = 0,title=0,removespike=True):  
    """
    This function plots either: 
    a)power spectra in csv files that are specfied in the list "paths".
    b)the nth power spectrum of a spectrogram saved in csv files specified in list "paths"
    (Can't mix power spectra files and dynamic spectra files in one list if nth_power_spec =/= 0)
    Assumes all spectra have same freq. range when setting axes.
    "norm":toggles normalizing all spectra with respect to peak PSD value
    "nth_power_spec":plots the nth row of the dataframe, applied to any files in list "paths" that contain dynamic spectrum
    peaks: toggles peak finding
    removespike: toggles removing DC spike and side artefacts
    """
    #setting figure size
    f = plt.figure()
    f.set_figwidth(16)
    f.set_figheight(9)
    
    #getting data from file and plotting spectrum
    for file in paths:
        
        data = make_data_array(file, nth_power_spec)
        
        power = data[0] #PSD values
        f_lower = data[1]
        f_upper = data[2]
        freq = np.linspace(f_lower,f_upper,len(power))
        
        #removing dc offset and side artefacts:
        if removespike == True:
            power,freq = remove_dcoffset(power,freq,func="cubic")
            power = truncate(power)
            freq = truncate(freq)
            
        #normalisation:
        if norm == True:
            max_val = abs(max(power))
            power = np.array(power)/max_val
            
        #peak finder function
        if peakfinder == True:
            peak_freqs,peak_values = peak_finder(power,freq)
            #plotting peaks found in peak finder function:
            plot_peaks(peak_freqs, peak_values)
        
        #plotting spectrum:
        plt.plot(freq,power, alpha = 0.7)
   
    #setting x ticks and labels:
    n_x = 9
    xticks_array = np.linspace(f_lower,f_upper,n_x)
    
    #rounding values for labels:
    xlabels = []
    for i in xticks_array:
        xlabels.append(np.round(i/1e6,2))
     
    xlabels = xlabels*u.MHz
    
    #Plotting
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("PSD (arbitrary units)")
    plt.xticks(xticks_array,xlabels)
    plt.legend(paths,loc = "best")
    #plt.legend(["HackRF One","LimeSDR mini"],loc = "best")
    if title == 0:
        title = paths[0]
    plt.title(title)
    plt.style.use("default")
    plt.grid()
    plt.show()

#argument parsing:
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("power_spec", help="name of file containing power spec to be plotted",type=str)
    parser.add_argument("-f2","--power_spec2", help="name of second file with power spec to be plotted, optional")
    parser.add_argument("-f3","--power_spec3",help="name of third file with power spec to be plotted, optional")
    parser.add_argument("--peaks",help="disable peak finder",action="store_true")
    parser.add_argument("-R","--spike_remover",help="disable dc spike remover",action="store_true")
    args = parser.parse_args()
    #optional arguments:
    #toggle peak finder:
    peak_bool = True
    if args.peaks:
        peak_bool = False
        
    print("peaks:",peak_bool)

    #toggle dc spike remover
    spike_bool = True
    if args.spike_remover:
        spike_bool = False

    print("spikes:",spike_bool)
        
    #plotting multiple spectra:
    paths = [args.power_spec]
    if args.power_spec2:
        paths.append(args.power_spec2)
    if args.power_spec3:
        paths.append(args.power_spec3)
    plot_power_spec(paths,peakfinder=peak_bool,removespike=spike_bool)
    
    

