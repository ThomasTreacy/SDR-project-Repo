#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:23:35 2022

@author: treacyth
"""

import numpy as np
import os 
import matplotlib.pyplot as plt
from astropy import units as u

os.chdir("D:\\python_files_usb")

import Soapy_power_plotter as spt
import Dynamic_spectra_plotter as dsp

filename = "D:\\SDR_backup\\SDR_data_files\\Week7\\limegrc1.bin"
data = np.fromfile(filename,dtype = np.float32)

print("\n",filename)
print("PSD values in file:",len(data))

#grc data files do not have nessecary parameters 
#saved in them, setting these manually below:
f_upper = 100e6
f_lower = 80e6
fft_bins = 2048
grc_av = 0 #number of averages performed in GRC

pow_spec = len(data)/fft_bins
print("pow_spec:",pow_spec)


i = 300# power spectrum to start from in data file
pyav = 320 #averages to perform in this script

lower = i*fft_bins
upper = (i+pyav)*fft_bins

data = data[lower:upper]
print("s:",data.shape)
print("data pow spec:",len(data)/fft_bins)


#shape into rows and cols
arr = np.reshape(data, (pyav,fft_bins))
print("shape:",arr.shape)

data = np.mean(arr, axis = 0)
print(data)
print(data.shape)

#Initialising plot:
fig = plt.figure()
fig.set_figwidth(14.4)
fig.set_figheight(8.1)

freq = np.linspace(f_lower,f_upper,fft_bins)
print(len(freq))

#Removing DC offset:
#data, freq = spt.remove_dcoffset(data,freq)

#Peak finding:
peak_freqs, peak_values = spt.peak_finder(data,freq,15)

#plotting peaks found in peak finder function:
spt.plot_peaks(peak_freqs, peak_values)

#freq axis ticks and labels:
x_ticks = np.linspace(f_lower,f_upper,9)
x_labels = [np.round(i/1e6)*u.MHz for i in x_ticks]

freq = spt.truncate(freq)
data = spt.truncate(data)

#plotting:
plt.plot(freq,data, alpha = 0.7)
t = "{} - 90 Mhz BW 20Mhz 2048bins GRC averages {}, Python averages {} window blackman".format(filename,grc_av,pyav)
plt.title(t + "\n")
plt.title("{} py averages {}".format(filename,pyav))

plt.ylabel("PSD (arbitrary units)")
plt.xlabel("Frequency (Mhz)")
plt.xticks(x_ticks,x_labels)
plt.grid()

#plt.legend(["soapy_power","GNU Radio companion"])
plt.show()


#hackgrc wk 8 90 Mhz BW 20Mhz 2048bins av 1000 from grc window blackman"