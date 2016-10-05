"""HW1 - DATA PROCESSING
Logan Halstrom
MAE 298 AEROACOUSTICS
HOMEWORK 1 - SIGNAL PROCESSING
CREATED: 04 OCT 2016
MODIFIY: 04 OCT 2016

DESCRIPTION: Read sound file of sonic boom and convert signal:
Narrow-band in Pa
Single-side power spectral density (FFT)
1/3 octave and octave band
"""

#IMPORT GLOBAL VARIABLES
from hw1_98_globalVars import *

import numpy as np
import pandas as pd

import os


import matplotlib.pyplot as plt



def ReadWav(filename):
    """Read a .wav file and return sampling frequency
    """
    from scipy.io import wavfile
    #Returns sample frequency and sampled data
    sampFreq, snd = wavfile.read(filename)
    return sampFreq, snd

#def SinglePowerSpec(fs, data):
#    """Find the single-sided power spectral density of a time function.
#    Requires an even amount of data points
#    fs --> sampling (source) frequency
#    data --> sampled data
#    """

#    #TIME
#    #calculate time of each signal, in seconds, from source frequency
#    N = len(data) #Number of data points in signal
#    dt = 1 / fs #time step
#    T = N * dt #total time interval of signal (s)
#    time = np.arange(N) * dt #individual sample times

#    #POWER SPECTRUM
#    fft = np.fft.fft(data) * dt #Fast-Fourier Transform
#    Sxx = np.abs(fft) ** 2 / T #Two-sided power spectrum

#    idx = range(int(N/2)) #Indices of single-sided power spectrum
#    Gxx = Sxx[idx] #Single-sided power spectrum
#    #Gxx = 2 * Sxx #Single-sided power spectrum

#    freqs = np.fft.fftfreq(data.size, dt) #Frequencies
#    #freqs = np.arange(N) / T #Frequencies

#    freqs = freqs[idx] #single-sided frequencies
#    print(freqs[:10])

def SPLi(P, Pref=20e-6):
    """Sound Pressure Level (SPL) in dB of a single pressue source (i)
    P --> pressure signal
    Pref --> reference pressure
    """
    PrmsSq = 0.5 * P ** 2 #RMS pressure squared
    return 10 * np.log10(PrmsSq / Pref ** 2)

def main(source):
    """input description

    """

    ####################################################################
    ### READ SOUND FILE ################################################
    ####################################################################

    df = pd.DataFrame() #Stores signal data
    #Read source frequency (fs) and signal in volts
    fs, df['V'] = ReadWav( '{}/{}'.format(datadir, source) )
    #Convert to pascals
    df['Pa'] = df['V'] * volt2pasc

    ####################################################################
    ### POWER SPECTRAL DENSITY #########################################
    ####################################################################

    #TIME
    #calculate time of each signal, in seconds, from source frequency
    N = len(df['Pa']) #Number of data points in signal
    dt = 1 / fs #time step
    T = N * dt #total time interval of signal (s)
    df['time'] = np.arange(N) * dt #individual sample times
    idx = range(int(N/2)) #Indices of single-sided power spectrum (first half)

    #POWER SPECTRUM
    fft = np.fft.fft(df['Pa']) * dt #Fast-Fourier Transform
    Sxx = np.abs(fft) ** 2 / T #Two-sided power spectrum
    #Gxx = Sxx[idx] #Single-sided power spectrum
    Gxx = 2 * Sxx[idx] #Single-sided power spectrum

    #FREQUENCY
    freqs = np.fft.fftfreq(df['Pa'].size, dt) #Frequencies
    #freqs = np.arange(N) / T #Frequencies
    freqs = freqs[idx] #single-sided frequencies

    #COMBINE POWER SPECTRUM DATA INTO DATAFRAME
    powspec = pd.DataFrame({'freq' : freqs, 'Gxx' : Gxx})




    #plt.figure()
    #plt.plot(df['time'], df['Pa'])
    #plt.show()








    ##NUMBER OF DISCRETE DATA POINTS
    #time = np.array(df['time'])
    #N = len(time)
    ##DATA TIME INTERVAL
    #TT = time[-1] - time[0]
    ##DATA TIME STEP
    #DT = time[1] - time[0]
    ##DATA BOUNDS
    #xmax = time[-1]
    #xmin = time[0]
    ##FFT OF DATA
    #fftfull = np.fft.fft(df['Pa'])
    ##Only use postitive frequencies (first half)
    #fft = fftfull[0:N/2-1]
    ##POWER SPECTRUM OF FFT
    #P = np.power(np.abs(fft)/(N/2.), 2)

    ##FREQUENCIES
    ##frequencies in hertz
    #freqh = np.fft.fftfreq(N, DT)
    #freqh = freqh[0:N/2-1]
    ##frequency numbers
    #freq = map(int, np.round(np.multiply(freqh, TT)))

    ##PLOT POWER SPECTRUM
    #plt.figure()
    #plt.plot(freqh, P)
    ###semi-log scale
    ##plt.set_yscale('log')
    ### ax.set_xscale('log')
    ### plt.xlim([0,0.25])
    ##plt.xlim([0,50])
    #plt.show()


    ##Power spectrum
    #fft = np.fft.fft(df['Pa']) * dt
    #Sxx = np.abs(fft) ** 2 / T
    #Gxx = 2 * Sxx

    #freqs = np.fft.fftfreq(df['Pa'].size, dt)
    ##freqs = np.arange(N) / T

    #print(freqs[:10])

    #idx = np.argsort(freqs)
    #idx = range(int(N/2))







    #plt.figure()
    #plt.plot(freqs, Gxx)
    #plt.xlim([0,50])
    #plt.show()

    ####################################################################
    ### FIND SOUND PRESSURE LEVEL IN dB ################################
    ####################################################################

    df['SPL'] = SPLi(df['Pa'])

    #plt.figure()
    #plt.plot(df['time'], df['SPL'])
    #plt.show()

    ####################################################################
    ### SONIC BOOM N-WAVE DURATION #####################################
    ####################################################################

    #Get times of maximum and minimum peaks of sonic boom
    tmax = float(df[df['Pa'] == max(df['Pa'])]['time']) #max press
    tmin = float(df[df['Pa'] == min(df['Pa'])]['time']) #min press
    #Duration of sonic boom N-wave (time from max pressure to min pressure)
    dt_Nwave = tmin - tmax


    ####################################################################
    ### SAVE DATA ######################################################
    ####################################################################

    #SAVE WAVE SIGNAL DATA
    df = df[['time', 'Pa', 'SPL', 'V']] #reorder
    df.to_csv( '{}/signal.dat'.format(datadir), sep=' ', index=False ) #save

    #SAVE POWER SPECTRUM DATA
    powspec.to_csv( '{}/powspec.dat'.format(datadir), sep=' ', index=False )

    #save single data (fs, Nwave)


if __name__ == "__main__":


    Source = 'Boom_F1B2_6.wav'

    main(Source)
