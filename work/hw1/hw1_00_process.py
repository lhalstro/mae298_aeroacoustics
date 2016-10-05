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

import numpy as np
import pandas as pd

import os


import matplotlib.pyplot as plt

#IMPORT GLOBAL VARIABLES
from hw1_98_globalVars import *

def ReadWav(filename):
    """Read a .wav file and return sampling frequency
    """
    from scipy.io import wavfile
    #Returns sample frequency and sampled data
    sampFreq, snd = wavfile.read(filename)
    return sampFreq, snd



def main(source):
    """input description

    """

    #READ SOUND FILE
    dat = pd.DataFrame() #Stores signal data
    #Read source frequency (fs) and signal in volts
    fs, dat['V'] = ReadWav( '{}/{}'.format(datadir, source) )
    T = 1 / fs #Period

    #calculate time of each signal, in seconds, from source frequency
    dat['time'] = np.arange( len(dat['V']) ) / fs

    #CONVERT TO PASCALS
    dat['Pa'] = dat['V'] * volt2pasc

    #plt.figure()
    #plt.plot(dat['time'], dat['Pa'])
    #plt.show()



    #Get times of maximum and minimum peaks of sonic boom
    tmax = float(dat[dat['Pa'] == max(dat['Pa'])]['time']) #max press
    tmin = float(dat[dat['Pa'] == min(dat['Pa'])]['time']) #min press
    #Duration of sonic boom N-wave (time from max pressure to min pressure)
    dt_Nwave = tmin - tmax




    #NUMBER OF DISCRETE DATA POINTS
    time = np.array(dat['time'])
    N = len(time)
    #DATA TIME INTERVAL
    TT = time[-1] - time[0]
    #DATA TIME STEP
    DT = time[1] - time[0]
    #DATA BOUNDS
    xmax = time[-1]
    xmin = time[0]
    #FFT OF DATA
    fftfull = np.fft.fft(dat['Pa'])
    #Only use postitive frequencies (first half)
    fft = fftfull[0:N/2-1]
    #POWER SPECTRUM OF FFT
    P = np.power(np.abs(fft) / (N/2.), 2.)

    #FREQUENCIES
    #frequencies in hertz
    freqh = np.fft.fftfreq(N, DT)
    freqh = freqh[0:N/2-1]
    #frequency numbers
    freq = map(int, np.round(np.multiply(freqh, TT)))

    #PLOT POWER SPECTRUM
    plt.figure()
    plt.plot(freq, P)
    #semi-log scale
    plt.set_yscale('log')
    # ax.set_xscale('log')
    # plt.xlim([0,0.25])
    plt.xlim([0,50])
    plt.show()




    ##Power spectrum
    #fft = np.fft.fft(dat['Pa'])
    #Sxx = np.abs(fft) ** 2 / T
    #Gxx = 2 * Sxx

    #freqs = np.fft.fftfreq(dat['Pa'].size, T)

    #print(freqs[:10])

    #idx = np.argsort(freqs)


    #plt.figure()
    #plt.plot(freqs[idx], Gxx[idx])
    #plt.show()



if __name__ == "__main__":


    Source = 'Boom_F1B2_6.wav'

    main(Source)
