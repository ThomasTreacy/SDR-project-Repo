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

def make_data_array(filename, file_format = "rtl_power"):
    """
    This function takes power spectrum in csv file "filename", extracts parameters and data
    Returns data as numpy array
    Also returns upper and lower freq. limits
    """
    if file_format == "rtl_power":
        #other file formats to be added
        spec_dataframe = pd.read_csv(filename,header = None)
        shape = spec_dataframe.shape
        
        #getting parameters 
        date = spec_dataframe.iloc[0,0]
        timestamp = spec_dataframe.iloc[0,1] 
        f_lower = spec_dataframe.iloc[0,2]
        f_upper = spec_dataframe.iloc[0,3]
        bin_size = spec_dataframe.iloc[0,4]
        samples = spec_dataframe.iloc[0,5]
        data = spec_dataframe.iloc[0,6:shape[1]]
    
    print("\n")
    print(filename)
    print("Timestamp:", date + timestamp)
    print("freq. range:",f_lower,f_upper)
    print("Bin size Hz:",bin_size)
    print("Samples:",samples)
    
    #storing data in numpy array
    data = np.array(data)
    return data, f_lower, f_upper
    

def plot_power_spec(paths,n_x,norm = False):
    """
    This function plots power spectra in csv files that are specfied in the list "paths".
    Assumes all spectra have same freq. range when setting axes.
    n_x is number of x axis ticks
    "norm" argument toggles normalizing all spectra with respect to peak PSD value
    false by default
    """
    #setting figure size
    f = plt.figure()
    f.set_figwidth(16)
    f.set_figheight(9)
    
    #getting data from file and plotting spectrum
    for file in paths:
        data = make_data_array(file)
        power = data[0] # array of PSD values arbitrary units
        f_lower = data[1] # array of freq values
        f_upper = data[2]
        
        freq = np.linspace(f_lower,f_upper,len(power))
        #normalisation:
        if norm:
            max_val = abs(max(power))
            power = np.array(power)/max_val
        
        plt.plot(freq,power)
        
    #setting x ticks and labels:
    xticks_array = np.linspace(f_lower,f_upper,n_x)
    
    #rounding values for labels:
    xlabels = []
    for i in xticks_array:
        xlabels.append(round(i/1e6,2))
     
    xlabels = xlabels*u.MHz
    
    #Plotting
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Power spectral density(arbitrary units)")
    plt.xticks(xticks_array,xlabels)
    plt.legend(paths)
    plt.title(paths[0])


#Plotting:
paths = ["hack_sig_gen1.csv","lime_sig_gen1.csv"]
plot_power_spec(paths,10)




            
        

