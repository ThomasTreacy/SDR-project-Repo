#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:23:58 2022

@author: treacyth

Originally created 13 10 22 to plot a spectrogram obtained using soapy_power.
This is a revised version started on 18 10 22
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.time import Time
from radiospectra.spectrogram2 import Spectrogram
from sunpy.net import attrs as a

import Soapy_power_plotter as spt
    
def load_spectrogram_data(file, file_format = "rtl_power",removespike=False):
    """
    Loads spectrogram data in the csv file "file". Imports file using pandas 
    and converts to numpy array. If multiple freq hops are in file these
    are organised such that one row of array has entire sweep of frequency spectrum
    outputs:
    data: spectrogram data
    freq: corresponding frequency array, y axis of spectrogram
    timestamps: time associated with each column of data, x axis of spectrogram
    """
    dyn_spec = np.array(pd.read_csv(file,header = None))
    
    #getting parameters and data from dataframe
    if file_format == "rtl_power":
        data = dyn_spec[0:,6:]
        bin_size = dyn_spec[0,4]
        print("raw data shape:",data.shape)
        
        print("last line:",dyn_spec[-1][:6])
        
        #organising separate freq hops:
            
        f_lowers = dyn_spec[0:,2]
        f_uppers = dyn_spec[0:,3]
        f_lower = min(f_lowers)
        f_upper = max(f_uppers)
        bins = data.shape[1]
        
        #get individual lower and upper freq values associated with hops:
        f_lowers_set = set()
        f_uppers_set = set()
        for i in range(len(f_lowers)):
            f_lowers_set.add(f_lowers[i])
            f_uppers_set.add(f_uppers[i])
        
        #converting sets to list so indexing works with them:
        f_lowers_set = list(f_lowers_set)
        f_uppers_set = list(f_uppers_set)
        
        #list of data arrays for each hop:
        hops_list = []
        timestamps = []
        #for timeseries:
        times = []
        
        #identifying a frequency hop by its lower freq. value
        for i1 in range(len(f_lowers_set)):
            hop_data = []#data for a given hop freq range
            for i in range(len(f_lowers)):
                #take all rows from data associated with given hop
                if f_lowers[i] == f_lowers_set[i1]:
                    #removing dc offset from row, need freq range:
                    freq = np.linspace(f_lowers_set[i1],f_uppers_set[i1],bins)
                    power = data[i]
                    if removespike == True:
                        power, freq = spt.remove_dcoffset(power,freq,func="linear")
                    hop_data.append(power)
                    #timestamps: take timestamps of lowest freq hop:
                    if f_lowers_set[i1] == f_lower:
                        date = dyn_spec[i,0]
                        time = dyn_spec[i,1]
                        timestamp = str(date) + str(time)
                        timestamps.append(timestamp)
                        times.append(time)
                
            hop_data = np.array(hop_data)
            hops_list.append(hop_data)
        
        #Join data from all hops in to one array for spectrogram:
        data = np.concatenate(hops_list, axis = 1)
        
        #putting dates and times into right format for astropy.time:
        timestamps = Time(timestamps)
        
        data = np.transpose(data)
    
        #when imported, all elements in array are of type "object",
        #converting data values to floats for plotting:
        data = data.astype(float)
        
    #print info:
    print("\n")
    print("filename:",file)    
    print("Bins:",bins)
    print("Freq.range:",f_lower,"Hz :",f_upper,"Hz")
    print("Bin size:",bin_size,"Hz")
    print("Sweeps:",data.shape[1])
    
    start_time = str(timestamps[0])[:19]
    end_time = str(timestamps[-1])[:19]
    print("Start time:",start_time)
    print("End time:  ",end_time)
    print("Average time for one sweep:",np.round(3585/data.shape[1],4),"seconds")
         
    #freq array:
    freq = np.linspace(f_lower,f_upper,data.shape[0])
    
    return data, freq, timestamps
        
def plot_timeseries(data, freq, timestamps, cent_freq):
    #Plot one timseries at a given freq f:
    #Inputs same as outputs of load_spectrogram_data, except cent_freq,
    #which is frequency at which timeseries is taken
    #Set relevant parameters:
    f_lower = freq[0]
    f_upper = freq[-1]
    bins = data.shape[0]
    start_time = str(timestamps[0])[:19]
    f = cent_freq
    n = int(bins*((f - f_lower)/(f_upper - f_lower)))#index of freq in array
    #take several rows from spectrogram around f and averge along f axis
    timeseries = data[n - 25 : n + 25,0:]
    print("ts_shape:",timeseries.shape)
    timeseries = np.mean(timeseries, axis = 0)
    
    #time axis ticks and labels:
    #(note: n_x and n_y may not be exact number of ticks 
    #due to floor division being used to calculate step)
    n_x = 10
    step = len(timestamps)//n_x
    x_ticks = np.arange(0,len(timestamps),step)
    x_labels = [str(timestamps[i])[10:19] for i in x_ticks]

    #plot the timeseries:
    plt.figure(figsize = [16,9],dpi = 500)
    plt.style.use("default")
    plt.plot(np.arange(len(timeseries)),timeseries)
    plt.xticks(x_ticks,x_labels, rotation = 60)
    plt.xlabel("Start time {} UTC".format(start_time))
    plt.ylabel("PSD (arbitrary units)")
    plt.title(file + " Time series at {} MHz".format(np.round(f/1e6,1)))
    #plt.ylim((-92,-82))
    plt.grid()
    plt.show()
    
def plot_pow_spec(data, freq, average=True, peakfinder=True,title = ""):
    """
    Plot one power spectrum from the spectrogram
    average: toggles averaging all of spectrogram into one power spectrum
    peakfinder: find peaks and plot onto power spectrum
    title: for figure.
    """
    #Set relevant parameters
    f_lower = freq[0]
    f_upper = freq[-1]   
    
    plt.figure(figsize = [16,9],dpi = 500)
    plt.style.use("default")
    #average enitre dynamic spectrum:
    if average == True:
        avs = data.shape[1]
        data = np.mean(data, axis = 1)
        plt.title(title  + " {} averages".format(avs))
    #take nth power spectrum:
    else:
        nth = 100
        data = data[0:,nth]
        plt.title(title + " power spectrum {}".format(str(nth)))
        
    #freq axis ticks and labels:
    x_ticks = np.linspace(f_lower,f_upper,9)
    x_labels = [np.round(i/1e6)*u.MHz for i in x_ticks]
        
    plt.plot(freq,data)
        
    #peakfinding:
    if peakfinder == True:
        peak_freqs, peak_vals = spt.peak_finder(data,freq)
        spt.plot_peaks(peak_freqs,peak_vals)
            
    plt.xticks(x_ticks,x_labels)
    plt.xlabel("Frequency (Mhz)")
    plt.ylabel("PSD (arbitrary units)")
    plt.grid()
    #plt.ylim((-127,-110))
    plt.show()

def plot_dyn_spec(data, freq, timestamps):
    """
    Plot spectrogram in the array "data", with frequency axis "freq" and time
    axis given by timestamps. Uses radiospectra.spectrogram2
    """
    freqs = freq*1e-6*u.MHz
    #Spectrogram2:
        
    #meta data
    meta = {
        'observatory': 'Monck',
        'instrument': 'LWA',
        'detector': 'LimeSDR mini',
        'freqs': freqs,
        'times': timestamps,
        'wavelength': a.Wavelength(freqs[0], freqs[-1]),
        'start_time': timestamps[0],
        'end_time': timestamps[-1]
    }
    spec = Spectrogram(data, meta)
    
    #plotting:
    plt.style.use("default")
    fig, ax = plt.subplots()
    plt.xlabel("Start time {} [UTC]".format(str(timestamps[0])[:19]))
    plt.ylabel("Frequency[MHz]")
    fig.set_figwidth(10)
    fig.set_figheight(7)
    spec.plot(ax,cmap = "viridis")#vmin = -100
    plt.show()

if __name__ == "__main__":
    file = "lime_lwa6.csv"
    data, freq, timestamps = load_spectrogram_data(file,removespike=False)
    plot_pow_spec(data, freq, title = file)
    plot_dyn_spec(data, freq, timestamps)
    plot_timeseries(data, freq, timestamps, 50e6)
    

    