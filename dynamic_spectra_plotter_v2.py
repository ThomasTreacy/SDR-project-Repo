#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:23:58 2022

@author: treacyth

Originally created 13 10 22 to plot dynamic spectra obtained using soapy_power.
This is a revised version created 18 10 22
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.time import Time
from radiospectra.spectrogram2 import Spectrogram
from sunpy.net import attrs as a
import soapy_power_test_v2 as spt

def plot_dynamic_spectrum(file, file_format = "rtl_power",plot_ts=False):
    """
    Plots a dynamic spectrum in the csv file "file". Imports file as
    pandas dataframe.
    Other arguments: 
    file_format: format of dynamic spectrum data in "file"
    plot_ts: toggles plotting a times series at a given freq within the spectrum
    """
    dyn_spec_dataframe = pd.read_csv(file,header = None)
    shape = dyn_spec_dataframe.shape # shape in form (rows, columns)
    print("shape:",shape)
    
    #getting parameters and data from dataframe
    if file_format == "rtl_power":
        #other file formats to be added
        data = dyn_spec_dataframe.iloc[0:shape[0],6:shape[1]]
        dates = np.array(dyn_spec_dataframe.iloc[0:shape[0],0])
        f_lower = dyn_spec_dataframe.iloc[0,2]
        f_upper = dyn_spec_dataframe.iloc[0,3]
        bin_size = dyn_spec_dataframe.iloc[0,4]
        times = np.array(dyn_spec_dataframe.iloc[0:shape[0],1])
        #putting dates and times into right format for astropy.time:
    
        timestamps = list(dates + times)
        timestamps = Time(timestamps)

        #putting axis right way round:
        data = data.T
        
    #freq array:
    freq = np.linspace(f_lower,f_upper,data.shape[0])
    
    cols = data.shape[1]
    
    #dc spike removal:
    list_of_series = []
    for i in range(cols):
        pow_spec = data.iloc[0:data.shape[0],i]
        pow_spec_arr = np.array(pow_spec)
        #print("l:",len(pow_spec),pow_spec_arr)
        pow_spec_arr,freq2 = spt.remove_dcoffset_and_sides(pow_spec_arr, freq, func = "linear")
        list_of_series.append(pow_spec_arr)
       
    freq = freq2
    data = pd.DataFrame(list_of_series)
    data = np.array(data)
    data = data.T
    #print("df:",data.shape)
    
    if plot_ts == True:
        #take one timseries at a given freq f:
        f = 85e6
        bins = int((f_upper - f_lower)/bin_size)
        n = int(bins*((f - f_lower)/(f_upper - f_lower)))
        timeseries = dyn_spec_dataframe.iloc[0:shape[0],n]
         
        #x axis:
        #(note: n_x and n_y may not be exact number of ticks 
        #due to floor division being used to calculate step)
        n_x = 10
        step = len(times)//n_x
        x_ticks = np.arange(0,len(times),step)
        x_labels = [times[i] for i in x_ticks]
    
        #plot the timeseries:
        plt.figure(figsize = [16,9],dpi = 500)
        plt.style.use("default")
        plt.plot(np.arange(len(timeseries)),timeseries)
        plt.xticks(x_ticks,x_labels, rotation = 60)
        plt.xlabel("Time UTC " + dates[0])
        plt.ylabel("PSD (arbitrary units")
        plt.title(file + " Time series at {} MHz".format(np.round(f/1e6,1)))
    else:
        #plotting dynamic spectrum:
        freqs = freq*1e-6*u.MHz
        #meta data
        meta = {
            'observatory': '',
            'instrument': 'Whip Antenna',
            'detector': 'LimeSDR mini',
            'freqs': freqs,
            'times': timestamps,
            'wavelength': a.Wavelength(freqs[0], freqs[-1]),
            'start_time': timestamps[0],
            'end_time': timestamps[-1]
        }
        plt.style.use("dark_background")
        spec = Spectrogram(data, meta )
        fig, ax = plt.subplots()
        fig.set_figwidth(10)
        fig.set_figheight(7)
        spec.plot(ax)

if __name__ == "__main__":
    #plotting entire dynamic spectrum:
    file = "lime_sun_spec3.csv"
    plot_dynamic_spectrum(file)
    #plotting the nth power spectrum from the dynamic spectrum:
    spt.plot_power_spec([file],nth_power_spec = 50,removespike=True,peakfinder=True)    
    
    