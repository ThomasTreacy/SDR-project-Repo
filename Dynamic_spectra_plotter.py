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
#from radiospectra.spectrogram2 import Spectrogram
from sunpy.net import attrs as a
import matplotlib.dates as mdates

os.chdir('D:\\python_files_usb')

import Soapy_power_plotter as spt

def bound(data, freq, timestamps, freq_bounds, time_bounds):
    """
    Limits the data to certain bounds for plotting in
    Inputs:
    data, freq, time same as output of load_spectrogram_data
    freq_bounds - list length two, lower then upper fraction of freq axis to bound
    time_bounds - list length two, lower then upper fraction of time axis to bound
    Each element of freq and time bounds lists must be between 0 and 1
    """
    #get indexes of bounds in arrays:
    freq_index1 = int((freq_bounds[0])*(data.shape[0]))
    freq_index2 = int((freq_bounds[1])*(data.shape[0]))
    time_index1 = int((time_bounds[0])*(data.shape[1]))
    time_index2 = int((time_bounds[1])*(data.shape[1]))
    #bound arrays using indexes:
    data = data[freq_index1:freq_index2,time_index1:time_index2]
    freq = freq[freq_index1:freq_index2]
    timestamps = timestamps[time_index1:time_index2]
    
    return data, freq, timestamps
    
def load_spectrogram_data(file, file_format = "rtl_power",removespike=False):
    """
    Loads spectrogram data in the csv file "file". Imports file using pandas 
    and converts to numpy array. If multiple freq hops are in file these
    are organised such that one row of array has entire sweep of frequency spectrum
    outputs:
    data: spectrogram data
    freq: corresponding frequency array, y axis of spectrogram
    timestamps: times associated with each column of data, astropy time objects
    """
    dyn_spec = np.array(pd.read_csv(file,header = None))
    
    #Several occasions where, with LimeSDR mini, one more hop at one freq range than another, correcting for this:
    #e.g. 14/12/2022 when observation was interupted half way through
    #if file == 'D:/SDRdata/Lime/From29Nov/12/08/2022_12_08_14%3A00.csv':
        #dyn_spec = dyn_spec[:-1]
    
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
        
        #get timestamp of first hop:
        timestamps = []
        #identifying a frequency hop by its lower freq. value
        for i1 in range(len(f_lowers_set)):
            for i in range(len(f_lowers)):
                if f_lowers[i] == f_lowers_set[i1]:
                    #timestamps: take timestamps of lowest freq hop:
                    if f_lowers_set[i1] == f_lower:
                        date = dyn_spec[i,0]
                        time = dyn_spec[i,1]
                        timestamp = str(date) + str(time)
                        timestamps.append(timestamp)
        
        #dc spike removal:  
        if removespike == True:
            data2 = []
            #freq array needed for spike remover func:
            freq = np.linspace(f_lower,f_upper,bins)
            for i in data:
                power, freq = spt.remove_dcoffset(i,freq, int_size = 10)
                data2.append(power)
            data2 = np.array(data2)
            print("data2:",data2.shape)
        
        #organise frequency hops
        hops = len(f_lowers_set)
        print("hops:",hops)
        
        hops_list = []
        
        for i in range(hops):
            hop_data = data[i::hops,:]
            hops_list.append(hop_data)
            print("hop:", hop_data.shape)
            
        #Join data from all hops in to one array for spectrogram:
        data = np.concatenate(hops_list, axis = 1)
        
        #putting dates and times into right format for astropy.time:
        timestamps = Time(timestamps)
        
        #Transpose so that freq is on y axis:
        data = np.transpose(data)
    
        #when imported, all elements in array are of type "object",converting data values to floats for plotting:
        data = data.astype(float)
        
        #freq array:
        freq = np.linspace(f_lower,f_upper,data.shape[0])
        
    #print info:
    print("\n")
    print("filename:",file)    
    print("Bins:",bins)
    print("Freq.range:",f_lower,"Hz :",f_upper,"Hz")
    print("Bin size:",bin_size,"Hz")
    print("Sweeps:",data.shape[1])
    
    return data, freq, timestamps
    
def plot_timeseries(data, freq, timestamps, cent_freq, top=None, SDR = "HackRF One"):
    #Plot one timseries at a given freq f:
    #Inputs same as outputs of load_spectrogram_data, except cent_freq,
    #which is frequency at which timeseries is taken
    print("Plotting timeseries...")
    #Set relevant parameters:
    f_lower = freq[0]
    f_upper = freq[-1]
    bins = data.shape[0]
    start_time = str(timestamps[0])[:19]
    f = cent_freq
    print(f)
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
    plt.xlabel("Start time {} [UTC]".format(start_time))
    plt.ylabel("PSD (arbitrary units)")
    plt.title("{} time series at {} MHz".format(SDR,np.round(f/1e6,1)))
    if top:
        plt.ylim(top = top)
    plt.grid()
    #plt.show()
    
def init_pow_spec(figsize = [13,7.34]):
    """
    initialises plot for power spectra - needed when plotting multiple power spectra in one graph
    """
    plt.figure(figsize = figsize,dpi = 500)
    plt.style.use("default")
    plt.grid()
    
def plot_pow_spec(data, freq, average=True,stdev =False,peakfinder=True,title = "",nth=100):
    """
    Plot one power spectrum from the spectrogram
    average: toggles averaging all of spectrogram into one power spectrum
    peakfinder: find peaks and plot onto power spectrum
    title: for figure.
    """
    print("Plotting power spectrum...")
    #Set relevant parameters
    f_lower = freq[0]
    f_upper = freq[-1]

    #average enitre dynamic spectrum:
    if average == True:
        avs = data.shape[1]
        data_for_plot = np.mean(data, axis = 1)
        plt.title(title  + " {} averages".format(avs))
    #take nth power spectrum:
    else:
        data_for_plot = data[0:,nth]
        plt.title(title + " power spectrum {}".format(str(nth)))
        
    #freq axis ticks and labels:
    x_ticks = np.linspace(f_lower,f_upper,9)
    x_labels = [np.round(i/1e6,2)*u.MHz for i in x_ticks]
  
    plt.plot(freq,data_for_plot)
    
    #plot plus and minus standard deviation curves:
    if stdev == True:
        stdev_data = np.std(data, axis = 1)
        plus_stdev = data_for_plot + stdev_data
        minus_stdev = data_for_plot - stdev_data
        plt.plot(freq,plus_stdev, alpha = 0.6)
        plt.plot(freq,minus_stdev, alpha = 0.6)
        plt.legend(["Average","Plus Std. Dev.","Minus Std. Dev."])
        
    #peakfinding:
    if peakfinder == True:
        peak_freqs, peak_vals = spt.peak_finder(data_for_plot,freq,n = 2)
        spt.plot_peaks(peak_freqs,peak_vals)
            
    plt.xticks(x_ticks,x_labels)
    plt.xlabel("Frequency (Mhz)")
    plt.ylabel("PSD (arbitrary units)")
    #plt.ylim((-127,-110))
    #plt.show()
    
    #callisto:
    """
    cal_data = np.array(pd.read_csv("D:/SDRData/birr_19dec/callisto.csv",header=None,delimiter=None))
    f_lower = 10e6
    f_upper = 115.5e6
    data = (cal_data[:,2] - 1050)/50 - 100
    freq = np.linspace(f_lower,f_upper,len(data))
    plt.plot(freq,data)
    """

def plot_dyn_spec(data, freq, timestamps, SDR, vmax = -60, vmin = -100):
    """
    Plot spectrogram in the array "data", with frequency axis "freq" and time
    axis given by timestamps. Uses radiospectra.spectrogram2
    SDR: name of software defined radio peripheral being used
    vmin, vmax: min and max PSD values in spectrogram
    """
    print("Plotting...")
    freqs = freq*1e-6*u.MHz
   
    #print info
    start_time = timestamps[0]
    end_time = timestamps[-1]
    duration = np.round((end_time - start_time).sec,1)
    start_time = str(start_time)[:19]
    print("Start time:",start_time)
    print("Duration:",str(duration) + " seconds")
    print("Average time for one sweep:",np.round(duration/data.shape[1],4),"seconds")
    
    #Lifted below from radiospectra2.spectrogram, got rid of data[:-1,:-1], which was causing me issues:
    fig, axes = plt.subplots(figsize = [9.72,7.29])
    axes.pcolormesh(timestamps.datetime, freqs.value, data, shading="auto", vmin = vmin, vmax = vmax)
    axes.set_xlim(timestamps.datetime[0], timestamps.datetime[-1])
    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    axes.xaxis.set_major_locator(locator)
    fig.autofmt_xdate()
    
    #Formatting:
    plt.xlabel("Start time {} [UTC]".format(str(timestamps[0])[:19]))
    plt.ylabel("Frequency[MHz]")
    #plt.title("RSTO,LWA,{}".format(SDR))
    ax = plt.gca()
    ax.invert_yaxis()

#%%
if __name__ == "__main__":
    #SDR = "HackRF One"
    SDR = "LimeSDR Mini"
    #os.chdir('D:/SDRdata/Lime/From29Nov/12/06')
    file = 'D:/SDRdata/hack/2022/12/30/2022_12_30_1500.csv'
    #file = 'D:/SDR_backup/SDR_data_files/Week6/lime_sat2.csv'
    data, freq, timestamps = load_spectrogram_data(file,removespike=True)
    #%%
    #data2, freq2, timestamps2 = load_spectrogram_data(file2,removespike=True)
    #Bounding freq and time:
    bgdata, bgfreq, bgtimestamps = bound(data, freq, timestamps,freq_bounds=[40/80,60/80],time_bounds=[(0/60),(3.50/60)])
    data2, freq2, timestamps2 = bound(data, freq, timestamps,freq_bounds=[40/80,60/80],time_bounds=[(26.50/60),(30/60)])
    
    data3 = data2 - bgdata
    
    plot_dyn_spec(data3, freq2, timestamps2, SDR=SDR, vmin = 0, vmax = 5)
    

    #%%
    nth = 10
    init_pow_spec()
    plot_pow_spec(data3, freq2, title = file,stdev=False,average=False,peakfinder=False,nth = nth)
    #plt.title("D:/SDRData/hack/2022/12/30/2022_12_30_1500.csv power spec 425 15:26:30.000")
    print("ps_time:", timestamps2[nth])
    
    #plot_timeseries(data2, freq2, timestamps2, cent_freq=38e6, SDR = SDR)
    #%%
    #plot_pow_spec(data, freq, title = file,stdev=False,average=True,peakfinder=False)
    plot_pow_spec(data2, freq2, title = file,stdev=False,average=True,peakfinder=False,nth = 3500)
    #plt.legend(["LimeSDR Mini","HackRF One"],loc = 'best')
    
#%%
    #plt.savefig(fname = "D:\Pics4report{}.png".format(), dpi = 200, bbox_inches = 'tight')



    