#!/bin/bash

#For regular LWA observations:

echo starting

#On Dell NUC
echo 'activating env'
source /home/monckobs/anaconda3/bin/activate soapy_env

#on laptop:
#source /home/treacyth/anaconda3/bin/activate sdr_project_env

var=$(date -u +%Y_%m_%d_%H:%M)

echo "$var" 

var2='debug_'
var4="$var2$var"

export PYTHONPATH=/usr/local/lib/python3.10/site-packages

#make folder for today's observations:
day=$(cut -c 9-10 <<< "$var")
month=$(cut -c 6-7 <<< "$var")
year=$(cut -c 1-4 <<< "$var")

directory=/home/monckobs/SDR_data_Nuc

#Check if there's folder for year, if not make it:
if [ -d "$directory/$year" ] 
then
    echo "folder exists"
else
    echo "folder doesn't exist, making it"
    mkdir "$directory/$year"
fi

#Check if there's folder for month, if not make it:
if [ -d "$directory/$year/$month" ] 
then
    echo "folder exists"
else
    echo "folder doesn't exist, making it"
    mkdir "$directory/$year/$month"
fi

#Check if there's folder for day, if not make it:
if [ -d "$directory/$year/$month/$day" ] 
then
    echo "folder exists"
else
    echo "folder doesn't exist, making it"
    mkdir "$directory/$year/$month/$day"
fi

directory="$directory/$year/$month/$day"

echo "$directory"
#HackRFOne:
soapy_power -e 3585 -f10M:90M -r20M -n2000 -b512 -G LNA=32,VGA=4,AMP=0 -D none --fft-window blackman -O $directory/$var.csv >& $directory/$var4.txt 

echo 'soapy_power -e 3585 -f10M:90M -r20M -n512 -b1024 -G LNA=32,VGA=4,AMP=0 -D none --fft-window blackman -O ' >> $directory/$var4.txt 

#Lime:
#soapy_power -e3585 -f15M:90M -r30M -n1000 -b1024 -G LNA=30,PGA=0,TIA=0 -D none --fft-window blackman -O $directory/$var.csv >& $directory/$var4.txt
