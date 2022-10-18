# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
import os, subprocess
import time
import sys
#from astropy.time import Time
from datetime import datetime
import pandas as pd
import io
import numpy as np

"""

#import pandas as pd

import subprocess
output = subprocess.run("soapy_power -r 20M -f 80M:100M -B 10k -O hack_thurs3.txt -T 1 -R -D none --fft-window boxcar",shell = True)
#subprocess.check_output
#subprocess.call
#subprocess.run
#print(output)

subprocess.run("dir", shell = True)

#%%

#1/10/22

#other modules that are needed:
import os
import datetime
import seaborn as sns
import matplotlib.pyplot as plt

my_df = pd.DataFrame({'A': [1, 2, 3], 'B': [4,5,6], 'C': [7,8,9], 'D': [10,11,12]})
my_df_transposed = my_df.T
print(my_df)
#print("transposed:","\n",my_df_transposed)
s = my_df.shape
print(s)
df2 = my_df.iloc[0:3,0:4]

print(df2)

df3 = my_df
df3.columns = ["H","I","J","K"]
print("\n",df3)

print("this should be 8:","\n",df3.iloc[1,2])

#%%
#importing using pandas test.csv
#need to do this with os first:
filename1=os.getcwd()+'/test.csv'
#first bit tells you directory of file - the current directory

#hackrf_df: pd.DataFrame = pd.read_csv(filename1, header=None, low_memory=False)
# Index sweeps by datetime
#datetime = pd.DatetimeIndex(hackrf_df[0])
#print(ArithmeticErrordatetime)#that error mes. appeared in the script!

#ax = sns.heatmap(data = hackrf_df.T) the csv file has headings and such that
#need to be sorted before plotting

ax = sns.heatmap(data = my_df,xticklabels = [1,2,3,4])
#plt.xticks([1,2,3])
plt.show()

#%%
print("my_df.index:",len(my_df.index))
print("len my df:",len(my_df))

#groupby func
def my_func(x):
    return(2*x)

#must give a func. or array
my_df = my_df.groupby(my_func)
print("test:",my_df)

    
#%%

print(os.getcwd())

#%%

f = open("test_text_doc.txt")
print(f.readline())
print(f.readline())


#%%
subprocess.run("date", shell = True)

#stop script for arg seconds
t.sleep(2)

subprocess.run("date", shell = True)

def function(x):
    return 2*x

#%%

class my_class:
    def my_func(self):
        print("test:",self.colour)
        print("test2:",self.house)
        """
        my_state = {
            'pwr_array': None
        }
        return my_state
    
    """
    #so that didn't really work
    def __init__(self,property1):
        self.colour = property1
        self.house = "car"

ob1 = my_class("blue")
ob1.my_func()    


#so that didn't really work

#%%

import numpy as np

test = np.array([1,2,3])
test2 = np.copy(test)
print(test2)e

#%%
#numpy zeros refresher
shape = np.array([10,10])
dtype = None
order = 'C'
z = np.zeros(shape,dtype,order)
print(z)

timestamps = [i for i in range(2000)]

n = 3
step = 2000//n
l = np.arange(0,2000,step)
lticks = l*(10/l[len(l) - 1])
print(l[0:20])
print("length:",len(l))

plt.xticks(lticks,l)
plt.imshow(z)




