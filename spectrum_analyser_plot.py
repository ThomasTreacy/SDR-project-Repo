# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 10:38:10 2022

@author: Thomas Treacy
#Loads csv files created with the Keysight N9340 software, can plot power spectrum or spectrogram with Dynamic_spectra_plotter.py

"""

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import Dynamic_spectra_plotter as dsp

def load_data(file):
    file_arr = np.array(pd.read_csv(file,header = None))
    
    num = int(file_arr[17,1])
    print("num:",num)
    
    #PSD values:
    data = file_arr[22:22 + num*(2):2,:-1]
    print("data:",data.shape)
    
    #timestamps:
    col4 = file_arr[21:21+ num*(2):2,3]
    print(col4)
    print("time:",col4.shape)
    
    #freq array:
    f_upper = float(file_arr[21,9])
    f_lower = float(file_arr[21,7])
    
    print(f_lower)
    print(type(f_lower))
    
    freq = np.linspace(f_lower,f_upper,data.shape[1])
    
    #transpose data so freq is on y axis:
    data = np.transpose(data)
    data = data.astype(float)
    return data, freq


data,freq = load_data('45_60_20 av.csv')
dsp.init_pow_spec()
dsp.plot_pow_spec(data,freq,average = True)
