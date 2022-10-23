#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:23:58 2022

@author: treacyth

Originally created 13 10 22 to plot dynamic spectra obtained using soapy_power.
This is a revised version created 18 10 22
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from astropy import units as u

import soapy_power_test_v2 as spt

def plot_dynamic_spectrum(file, n_x, n_y,file_format = "rtl_power"):
    """
    Plots a dynamic spectrum in the csv file "file". Imports file as
    pandas dataframe.
    Other arguments: 
    n_x: rough number of x axis ticks
    n_y: rough number of y axis ticks
    file_format: format of dynamic spectrum data in "file"
    """
    dyn_spec_dataframe = pd.read_csv(file,header = None)
    shape = dyn_spec_dataframe.shape # shape in form (rows, columns)
    
    #getting parameters and data from dataframe
    if file_format == "rtl_power":
        #other file formats to be added
        date = dyn_spec_dataframe.iloc[0,0]
        f_lower = dyn_spec_dataframe.iloc[0,2]
        f_upper = dyn_spec_dataframe.iloc[0,3]
        bin_size = dyn_spec_dataframe.iloc[0,4]
        timestamps = dyn_spec_dataframe.iloc[0:shape[0],1]
        #storing timestamps in numpy array:
        time_arr = np.array(timestamps)
        
        #isolating data in data frame
        data = dyn_spec_dataframe.iloc[0:shape[0],6:shape[1]]
        
        #putting axis right way round:
        data = data.T

    #x axis:
    
    #(note: n_x and n_y may not be exact number of ticks 
    #due to floor division being used to calculate step)
    step = len(timestamps)//n_x
    
    x_ticks = np.arange(0,len(timestamps),step)
    x_labels = [timestamps[i] for i in x_ticks]
    
    #y axis:
    freq = np.linspace(f_lower,f_upper,data.shape[0])
    
    step = len(freq)//n_y
    y_ticks = np.arange(0,data.shape[0],step)
    y_labels = [freq[i] for i in y_ticks]
    
    #rounding floats for labels and adding freq. MHz
    ylabels2 = []
    for i in y_labels:
        ylabels2.append(round(i*1e-6,1)*u.MHz)
    y_labels = ylabels2
    
    #plotting:
    plt.figure(figsize = [8,4.5],dpi = 500)
    ax = sns.heatmap(data, cbar_kws={'label': 'PSD (arbitrary units)'})
    plt.xticks(x_ticks,x_labels, rotation = 60)
    plt.yticks(y_ticks,y_labels)
    plt.title(file)
    plt.ylabel("Frequency (MHz)")
    plt.xlabel("Time UTC " + date)

#plotting entire dynamic spectrum:
file = "lime_sat_spec30min.csv"
plot_dynamic_spectrum(file,10,10)

#plotting the nth power spectrum from the dynamic spectrum:
spt.plot_power_spec([file],10,nth_power_spec = 0)    
