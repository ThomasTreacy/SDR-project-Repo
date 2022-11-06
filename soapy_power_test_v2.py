#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:25:33 2022

@author: treacyth

Originally created 4 10 22 to plot two power spectra from a txt or csv file obtained using soapy_power. 
This is a revision started on 18 10 22
Data files must be in current working directory, must be csv file

"""


import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def peak_finder(power,freq):
    """
    Uses scipy.signal.find_peaks to find peaks in power array
    freq: array of freq values corresponding to power array
    """
    peak_indices = find_peaks(power,prominence=2)[0]
    print("testing")
    #getting PSD values and freq of peaks:
    peak_freqs = []
    peak_values = []
    for i in peak_indices:
        peak_values.append(power[i])
        peak_freqs.append(freq[i])
    return np.array(peak_freqs),np.array(peak_values)

def make_data_array(filename, nth_power_spec = 0,file_format = "rtl_power"):
    """
    This function takes a csv file "filename" containing either power or dynamic spectrum, 
    extracts parameters and plots ONE power spectrum as numpy array
    Also returns upper and lower freq. limits
    "nth_power_spec": plot the nth power spectrum if filename contains dynamic spectrum
    """
    if file_format == "rtl_power":
        #other file formats to be added
        spec_dataframe = pd.read_csv(filename,header = None)
        shape = spec_dataframe.shape
        
        #getting parameters and data
        date = spec_dataframe.iloc[0,0]
        timestamp = spec_dataframe.iloc[0,1] 
        f_lower = spec_dataframe.iloc[0,2]
        f_upper = spec_dataframe.iloc[0,3]
        bin_size = spec_dataframe.iloc[0,4]
        samples = spec_dataframe.iloc[0,5]
        data = spec_dataframe.iloc[nth_power_spec,6:shape[1]]
        #data: power spectral density (PSD) values
    
    print("\n")
    print(filename)
    print("Timestamp:", date + timestamp)
    print("freq. range:",f_lower,f_upper)
    print("Bin size Hz:",bin_size)
    print("Samples:",samples)
    print("\n")
    
    #storing data in numpy array
    data = np.array(data)
    return data, f_lower, f_upper

def poly_order_3(x,c,d,e,f):
    """
    Cubic function for curve ftting
    """
    return c*x**3 + d*x**2 + e*x + f

def remove_dcoffset_and_sides(power,freq,func="cubic"):
    """
    This function returns array that is the input array with the central values replaced,
    to remove DC spike, and values at edge of input array cut off to remove side artefacts
    power: array of PSD values
    freq: array of freq values corresponding to "power"
    func: kind of fitting function used to replace DC spike, either cubic or linear
    """
    
    num = 20#size of interval
    #defining dc spike region:
    lower = int(len(power)/2 - num)
    upper = int(len(power)/2 + num)
    
    #getting data points either side of DC spike gap
    num = 15
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
    
    
    #removing side artefacts
    num = 15
    lower = num
    upper = int(len(power) - num)
    power = power[lower:upper] 
    freq = freq[lower:upper]
    
    return power, freq

def plot_power_spec(paths,peakfinder = True,norm = False,nth_power_spec = 0,title=0,removespike=True):  
    """
    This function plots either: 
    a)power spectra in csv files that are specfied in the list "paths".
    b)the nth power spectra of a dynamic spectrum saved in csv files specified in list "paths"
    (Can't mix power spectra files and dynamic spectra files in one list)
    Assumes all spectra have same freq. range when setting axes.
    "norm":toggles normalizing all spectra with respect to peak PSD value
    "nth_power_spec":plots the nth row of the dataframe, applied when files in list "paths" that contain dynamic spectum
    peakfinder: toggles peak finding
    title: if not specified, title of plot is first file in paths
    removespike: toggle removing dc spike and side artefacts
    """
    #setting figure size
    f = plt.figure()
    f.set_figwidth(16)
    f.set_figheight(9)
    
    #getting data from file and plotting spectrum
    for file in paths:
        data = make_data_array(file, nth_power_spec)
        power = data[0] #Power spec
        
        # freq array:
        f_lower = data[1] 
        f_upper = data[2]
        freq = np.linspace(f_lower,f_upper,len(power))
        
        #removing dc offset and side artefacts:
        if removespike == True:
            power,freq = remove_dcoffset_and_sides(power,freq,func="cubic")
        
        #normalisation:
        if norm == True:
            max_val = abs(max(power))
            power = np.array(power)/max_val
            
        #peak finder function
        if peakfinder == True:
            peaks = peak_finder(power,freq)
            #plotting peaks found in peak finder function:
            x = peaks[0]
            y = peaks[1]
            for i in range(len(x)):
                #freqs in x are in Hz - converting to Mhz for labelling
                x_label = str(np.round(x[i]/1e6,2)*u.MHz)
                plt.text(x[i],y[i],x_label + "\n" + str(y[i]))
        
        #plotting spectrum:
        plt.plot(freq,power, alpha = 0.7)
   
    #setting x ticks and labels:
    n_x = 9
    xticks_array = np.linspace(f_lower,f_upper,n_x)
    
    #rounding values for labels:
    xlabels = []
    for i in xticks_array:
        xlabels.append(np.round(i/1e6,2)*u.MHz)
    
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
    

if __name__ == "__main__":
    paths = ["lime_sat1.csv"]
    plot_power_spec(paths,peakfinder=True,removespike=True)


