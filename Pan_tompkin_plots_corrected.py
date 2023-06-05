# -*- coding: utf-8 -*-                                                        
"""
Created on Thu May 25 17:32:37 2023

@author: Asus
"""

WINDOWSIZE = 20  
NOSAMPLE = -32000
BUFFSIZE = 600
DELAY = 22  
FS=360     #SAMPLING INTERVAL, NOT FREQUENCY

import os

os.chdir("C:\\Users\\Asus\\Documents\\Python Programs")

import pandas as pd

'''def output(op,input_val):  
    if(op==True):
        print (input_val," - ",op)'''

signal = [0] * BUFFSIZE
dcblock = [0] * BUFFSIZE
lowpass = [0] * BUFFSIZE
highpass = [0] * BUFFSIZE
derivative = [0] * BUFFSIZE
squared = [0] * BUFFSIZE
integral = [0] * BUFFSIZE
outputSignal = [0] * BUFFSIZE

rr1 = [0] * 8
rr2 = [0] * 8
rravg1 = 0
rravg2 = 0
rrlow = 0
rrhigh = 0
rrmiss = 0

i = 0
j = 0
sample = 0
lastQRS = 0
lastSlope = 0
currentSlope = 0

current=0
peak_i=0
peak_f=0
threshold_i1=0
threshold_i2=0
threshold_f1 =  0
threshold_f2 = 0
spk_i = 0 
spk_f = 0
npk_i = 0 
npk_f = 0
qrs = False
regular = True
prevRegular = False

for i in range(8):
    rr1[i]=0
    rr2[i]=0

signal_s = [] 
dcblock_s = [] 
lowpass_s = [] 
highpass_s = []
derivative_s = [] 
squared_s = [] 
integral_s = []

###############################################################################
def average_until_zero(lst):
    total_sum = 0
    count = 0

    for num in lst:
        if num == 0:
            break
        total_sum += num
        count += 1

    if count == 0:
        return 0

    average = total_sum / count
    return average

import statistics

def stddev_until_zero(lst):
    non_zero_elements = []

    for num in lst:
        if num == 0:
            break
        non_zero_elements.append(num)

    if len(non_zero_elements) == 0:
        return 0  

    standard_deviation = statistics.stdev(non_zero_elements)
    return standard_deviation

import numpy as np
from scipy.signal import periodogram 

def psd_until_zero(lst):
    signal_array = np.array(signal)
    sampling_interval = FS  
    zero_index = np.where(signal_array == 0)[0]
    if len(zero_index) == 0:
        zero_index = len(signal_array)
    else:
        zero_index = zero_index[0]
    sampling_frequency = 1000 / sampling_interval
    frequency = np.fft.rfftfreq(zero_index, d=1/sampling_frequency)
    f, psd = periodogram(signal_array[:zero_index], fs=sampling_frequency)
    return frequency, psd

from scipy import stats

def skewness_until_zero(data):
    
    zero_index = np.where(data == 0)[0]
    if len(zero_index) == 0:
        zero_index = len(data)
    else:
        zero_index = zero_index[0]
    
    skewness = stats.skew(data[:zero_index])
    return skewness

def kurtosis_until_zero(data):
    # Find the index where the data becomes zero
    zero_index = np.where(data == 0)[0]
    if len(zero_index) == 0:
        zero_index = len(data)
    else:
        zero_index = zero_index[0]
    
    # Calculate the kurtosis until the zero index
    kurtosis = stats.kurtosis(data[:zero_index])
    return kurtosis

###############################################################################
    
def PanTompkins(input_val):
    
    global rr1 
    global rr2 
    global rravg1 
    global rravg2 
    global rrlow 
    global rrhigh 
    global rrmiss 

    global i
    global j 
    global sample 
    global lastQRS 
    global lastSlope 
    global currentSlope 

    global current
    global peak_i
    global peak_f
    global threshold_i1
    global threshold_i2
    global threshold_f1  
    global threshold_f2 
    global spk_i  
    global spk_f 
    global npk_i  
    global npk_f 
    global qrs 
    global regular 
    global prevRegular
    
    global signal_s
    
    if sample>=BUFFSIZE:
        for i in range(BUFFSIZE - 1):
            signal[i] = signal[i + 1]
            dcblock[i] = dcblock[i + 1]
            lowpass[i] = lowpass[i + 1]
            highpass[i] = highpass[i + 1]
            derivative[i] = derivative[i + 1]
            squared[i] = squared[i + 1]
            integral[i] = integral[i + 1]
            outputSignal[i] = outputSignal[i + 1]
        current = BUFFSIZE - 1   
    else:
        current=sample
    

    
    signal[current] = input_val
    signal_s.append(signal[current])
    
    
    sample+=1
    
    if current >= 1:
        dcblock[current] = signal[current] - signal[current-1] + 0.995*dcblock[current-1]
    else:
        dcblock[current] = 0
    
    dcblock_s.append(dcblock[current])
    lowpass[current] = dcblock[current]
    
    if current >= 1:
        lowpass[current] += 2 * lowpass[current-1]
    if current >= 2:
        lowpass[current] -= lowpass[current-2]
    if current >= 6:
        lowpass[current] -= 2 * dcblock[current-6]
    if current >= 12:
        lowpass[current] += dcblock[current-12]
    
    lowpass_s.append(lowpass[current])
    
    highpass[current] = -lowpass[current]
    
    if current >= 1:
        highpass[current] -= highpass[current-1]
    if current >= 16:
        highpass[current] += 32 * lowpass[current-16]
    if current >= 32:
        highpass[current] += lowpass[current-32]
    
    highpass_s.append(highpass[current])
    
    derivative[current] = highpass[current]
    
    if current > 0:
        derivative[current] -= highpass[current-1]
    
    derivative_s.append(derivative[current])
    
    squared[current] = derivative[current] * derivative[current]
    squared_s.append(squared[current])
    
    integral[current] = 0
    for i in range(WINDOWSIZE):
        if current >= i:
            integral[current] += squared[current - i]
        else:
            break
    
    
    integral[current] /= i
    qrs = False
    
    integral_s.append(integral[current])
    
    
    if integral[current] >= threshold_i1 or highpass[current] >= threshold_f1:
        peak_i = integral[current]
        peak_f = highpass[current]
    
    if integral[current] >= threshold_i1 and highpass[current] >= threshold_f1:
        if sample > lastQRS + FS/5:
            if sample <= lastQRS + int(0.36*FS):
                currentSlope = 0
                for j in range(current - 10, current + 1):
                    if squared[j] > currentSlope:
                        currentSlope = squared[j]
    
                if currentSlope <= lastSlope/2:
                    qrs = False
                else:
                    spk_i = 0.125 * peak_i + 0.875 * spk_i
                    threshold_i1 = npk_i + 0.25 * (spk_i - npk_i)
                    threshold_i2 = 0.5 * threshold_i1
    
                    spk_f = 0.125 * peak_f + 0.875 * spk_f
                    threshold_f1 = npk_f + 0.25 * (spk_f - npk_f)
                    threshold_f2 = 0.5 * threshold_f1
    
                    lastSlope = currentSlope
                    qrs = True
            else:
                currentSlope = 0
                for j in range(current - 10, current + 1):
                    if squared[j] > currentSlope:
                        currentSlope = squared[j]
                spk_i = 0.125 * peak_i + 0.875 * spk_i
                threshold_i1 = npk_i + 0.25 * (spk_i - npk_i)
                threshold_i2 = 0.5 * threshold_i1
    
                spk_f = 0.125 * peak_f + 0.875 * spk_f
                threshold_f1 = npk_f + 0.25 * (spk_f - npk_f)
                threshold_f2 = 0.5 * threshold_f1

                
                lastSlope = currentSlope
                qrs = True
                
        else:
            peak_i = integral[current]
            npk_i = 0.125*peak_i + 0.875*npk_i
            threshold_i1 = npk_i + 0.25*(spk_i - npk_i)
            threshold_i2 = 0.5*threshold_i1
            peak_f = highpass[current]
            npk_f = 0.125*peak_f + 0.875*npk_f
            threshold_f1 = npk_f + 0.25*(spk_f - npk_f)
            threshold_f2 = 0.5*threshold_f1
            qrs = False
            outputSignal[current] = qrs
            #if sample > DELAY + BUFFSIZE:
            #    output(outputSignal[0],)
            
    if qrs:
        rravg1=0
        for i in range(7):
            rr2[i]=rr1[i+1]
            rravg1+rr1[i]
            
        rr1[7] = sample - lastQRS
        lastQRS = sample
        rravg1 += rr1[7]
        rravg1 *= 0.125
        
        if rr1[7] >= rrlow and rr1[7] <= rrhigh:
            rravg2 = 0
            for i in range(7):
                rr2[i] = rr2[i+1]
                rravg2 += rr2[i]
            rr2[7] = rr1[7]
            rravg2 += rr2[7]
            rravg2 *= 0.125
            rrlow = 0.92 * rravg2
            rrhigh = 1.16 * rravg2
            rrmiss = 1.66 * rravg2
            
        prevRegular = regular
        if rravg1 == rravg2:
            regular = True
        else:
            regular = False
            if prevRegular:
                threshold_i1 /= 2
                threshold_f1 /= 2
    
    else:
        if sample - lastQRS > rrmiss and sample > lastQRS + FS/5:
            for i in range(current - (sample - lastQRS) + FS//5, current):

                if (integral[i] > threshold_i2) and (highpass[i] > threshold_f2):
                    currentSlope = 0
                    for j in range(i - 10, i + 1):
                        if squared[j] > currentSlope:
                            currentSlope = squared[j]
                    if currentSlope < lastSlope/2 and i + sample < lastQRS + 0.36*lastQRS:
                        qrs = False
                    else:
                        peak_i = integral[i]
                        peak_f = highpass[i]
                        spk_i = 0.25*peak_i + 0.75*spk_i
                        spk_f = 0.25*peak_f + 0.75*spk_f
                        threshold_i1 = npk_i + 0.25*(spk_i - npk_i)
                        threshold_i2 = 0.5*threshold_i1
                        lastSlope = currentSlope
                        threshold_f1 = npk_f + 0.25*(spk_f - npk_f)
                        threshold_f2 = 0.5*threshold_f1
                        
                        # RR Average 1
                        rravg1 = 0
                        for j in range(7):
                            rr1[j] = rr1[j+1]
                            rravg1 += rr1[j]
                        rr1[7] = sample - (current - i) - lastQRS
                        qrs = True
                        lastQRS = sample - (current - i)
                        rravg1 += rr1[7]
                        rravg1 *= 0.125
                        
                        if (rr1[7] >= rrlow) and (rr1[7] <= rrhigh):
                            rravg2 = 0
                            for i in range(7):
                                rr2[i] = rr2[i+1]
                                rravg2 += rr2[i]
                            rr2[7] = rr1[7]
                            rravg2 += rr2[7]
                            rravg2 *= 0.125
                            rrlow = 0.92 * rravg2
                            rrhigh = 1.16 * rravg2
                            rrmiss = 1.66 * rravg2
                    
                        prevRegular = regular
                        if rravg1 == rravg2:
                            regular = True
                        else:
                            regular = False
                            if prevRegular:
                                threshold_i1 /= 2
                                threshold_f1 /= 2
                        break
            if qrs:
                outputSignal[current] = False
                outputSignal[i] = True
                #if sample > DELAY + BUFFSIZE:
                #    output(outputSignal[0])
        
        if not qrs:
            if (integral[current] >= threshold_i1) or (highpass[current] >= threshold_f1):
                peak_i = integral[current]
                npk_i = 0.125*peak_i + 0.875*npk_i
                threshold_i1 = npk_i + 0.25*(spk_i - npk_i)
                threshold_i2 = 0.5*threshold_i1
                peak_f = highpass[current]
                npk_f = 0.125*peak_f + 0.875*npk_f
                threshold_f1 = npk_f + 0.25*(spk_f - npk_f)
                threshold_f2 = 0.5*threshold_f1
    
    outputSignal[current] = qrs
    #if sample > DELAY + BUFFSIZE:
        #output(outputSignal[0])
    
    #print("intgrl: ", integral[0]," i1: ",threshold_i1, " f1: ",threshold_f1)
    return qrs;
        

        
        
    #for i in range(1, BUFFSIZE):
    #    output(outputSignal[i])
    
    
data=pd.read_csv('sample.csv',header=None)
column=data.iloc[:,0]

count=0
truecount=0
peaksample=[]
for input_val in column:
    result=PanTompkins(int(input_val)) 
    
    if(result):
        truecount+=1
        #print(count," - ",input_val," - ",result)
        peaksample.append(count)
        
        
    count+=1
print ('Peak Samples using Pan Tompkins',peaksample)
    
mat=[]
MATRIXWIDTH=356
temp=[0]*MATRIXWIDTH
for i in range(truecount-1):
    mat.append(temp.copy())
    
    
curr_ind=0
next_ind=1
samplenumber=1
row=0
col=0
firstpeakfound=False
for input_val in column:
    if curr_ind < len(peaksample) and next_ind < len(peaksample):
        if(samplenumber>=peaksample[curr_ind] and samplenumber<peaksample[next_ind]):
            firstpeakfound=True
            mat[row][col]=int(input_val)
            col+=1
        else:
            if firstpeakfound==True and row<truecount-2:
                row+=1
                col=0
                mat[row][col]=int(input_val)
                col+=1
                curr_ind+=1
                next_ind+=1
                
        samplenumber+=1
        
              

#print(mat)
import matplotlib.pyplot as plt

from scipy.signal import find_peaks

peaks=find_peaks(integral_s,distance=40)
print('Using find_peaks() with distance=40',peaks)
peaks=find_peaks(integral_s)
#print('Using find_peaks() without any other parameter',peaks)              
plt.figure()
plt.plot(signal_s)
plt.xlabel('Sample Number')
plt.ylabel('Signal o/p')
plt.title('SIGNAL PLOT')

plt.figure()
plt.plot(dcblock_s)
plt.xlabel('Sample Number')
plt.ylabel('DC block o/p')
plt.title('DC BLOCK PLOT')

plt.figure()
plt.plot(lowpass_s)
plt.xlabel('Sample Number')
plt.ylabel('Lowpass o/p')
plt.title('LOWPASS PLOT')

plt.figure()
plt.plot(highpass_s)
plt.xlabel('Sample Number')
plt.ylabel('Highpass o/p')
plt.title('HIGHPASS PLOT')
                       
plt.figure()
plt.plot(derivative_s)
plt.xlabel('Sample Number')
plt.ylabel('Derivative o/p')
plt.title('DERIVATIVE PLOT')

plt.figure()
plt.plot(squared_s)
plt.xlabel('Sample Number')
plt.ylabel('Squared o/p')
plt.title('SQUARED PLOT')
                                     
plt.figure()
plt.plot(integral_s)
plt.xlabel('Sample Number')
plt.ylabel('Integral o/p')
plt.title('INTEGRAL PLOT')

plt.show()

i=0;
while(i<truecount-1):
    row=mat[i]
    avg=average_until_zero(mat[i])
    stddev=stddev_until_zero(mat[i])
    ind,psd=psd_until_zero(mat[i])
    skw=skewness_until_zero(mat[i])
    kurt=kurtosis_until_zero(mat[i])
    
    print("For Row ",i," ")
    print("    Mean: ",avg,"  Std Devation:  ",stddev)
    print("    Skewness:  ",skw,"Kurtosis:  ",kurt)
    print("    Freqeuncy: ",ind,"  PSD",psd)
    i+=1
   
    
