#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 17:19:19 2022

@author: treacyth

A script originally used to combine power spectra from soapy power into dynamic spectra and plot them. Slow because of the
time it took to repeatly stop and start soapy power. Became obsolete once it was realised that soapy power includes functionality for 
collecting power spectra.
"""

import subprocess
import matplotlib.pyplot as plt
import time as t
import numpy as np
from astropy import units as u
import os
import seaborn as sns
import SDR_project_funcs as sdr

#name for dynamic spectrum:
spec_name = "hack_nov_spec"
#making a folder in cwd for power spectra files:
first_dir = os.getcwd()
print(first_dir)

new_dir = spec_name + "_folder"
os.mkdir(new_dir)
os.chdir(first_dir + "/" + new_dir)

print(os.getcwd())
#changes back to first dir. at end of script

#number of power spectra obtained:
power_spec_num = 20
spec_files = []

#set freq. range in MHz:
f_lower = 80
f_upper = 100

sample_rate = f_upper - f_lower

#bin size kHz
bin_size = 10
    
#looping to get multiple power spectra:
for i in range(power_spec_num):
    filename = "{}_{}.txt"
    filename = filename.format(spec_name,"%s") % i
    spec_files.append(filename)
    
    string = "soapy_power -r {}M -f {}M:{}M -B {}k -O %s -T 1 -R -D none --fft-window boxcar" % filename
    string = string.format(sample_rate,f_lower,f_upper,bin_size)
    output = subprocess.run(string, shell = True)
    print("loop no:",i)
    print("string:",string)
 
#---------------------------------------------

os.chdir("/home/treacyth/Documents/LimeSDR_mini_sunday_folder")
#getting data from files for plotting:
#formatting assumes one hop per file, certain freq. ranges produce two hops in a file, need to find out why
array_list = []
time_stamps = []
 
 def make_data_list(filename):
    #this function takes text file and puts power values into list
    #and also makes list of corresponding freq. values
    
    #opening file:
    file = open(filename)
    contents = file.read()
    
    #converting string to list:
        
    contents_list = contents.split()
    date = contents_list[0]
    
    #print("test:", contents_list)
    
    #taking upper and lower freq. values stored in file:
    f_lower = float(contents_list[2].replace(",",""))
    f_upper = float(contents_list[3].replace(",",""))
    
    #getting actual data values of file:
    data_list = contents_list[6:]
    
    data_list_2 = []
    #removing commas from strings and making strings into floats:
    for i in data_list:
        data_list_2.append(float(i.replace(",", "")))
    
    #print("test:",data_list_2)
    
    #print("test:",f_lower,f_upper,len(data_list_2))
    
    #freq. list:
    n = np.linspace(f_lower,f_upper,len(data_list_2))
    
    print("\n")
    print(filename)
    print("Timestamp:",contents_list[0],contents_list[1])
    print("freq. range:",f_lower,f_upper)
    print("Bin size Hz:",contents_list[4])
    
    
    return data_list_2,n
    
for i in spec_files:
    m = make_data_list(i)
    powers = m[0]
    time = m[2]
    time_stamps.append(time)
    array_list.append(powers)
    
    
array = np.array(array_list)
#putting time and freqs. on correct axes:
array = array.T

#x axis labelling:
#labelling every second time stamp:
xticks_array = np.arange(0,len(time_stamps),2)

time_stamps2 = []
for i in range(0,len(time_stamps),2):
    time_stamps2.append(time_stamps[i])
time_stamps = time_stamps2
#print("test:",time_stamps)  

#y axis labelling
bins = len(array_list[0])
print("bins:",bins)
yticks_array = np.linspace(0,bins,10)
ylabels = np.linspace(f_lower,f_upper,10) 

#rounding floats
ylabels2 = []

for i in ylabels:
    ylabels2.append(round(i,1))
ylabels = ylabels2*u.MHz

plot = sns.heatmap(array)
plt.title(spec_name)
plt.xticks(xticks_array,time_stamps,rotation = 60)
plt.yticks(yticks_array,ylabels)
plt.xlabel("Time")
plt.ylabel("Frequency")
plt.show()

#%%
#changing back to initial dir:
os.chdir(first_dir)
print(os.getcwd())
    
    
    
