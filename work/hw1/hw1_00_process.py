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
    as a function of time
    P --> pressure signal
    Pref --> reference pressure
    """
    PrmsSq = 0.5 * P ** 2 #RMS pressure squared
    return 10 * np.log10(PrmsSq / Pref ** 2)

def SPLf(Gxx, T, Pref=20e-6):
    """Sound Pressure Level (SPL) in dB of a single pressue source (i)
    P --> pressure signal
    Pref --> reference pressure
    """
    return 10 * np.log10( (Gxx / T) / Pref ** 2 )

def OctaveBounds(fc, octv=1):
    """Get upper/lower frequency bounds for given octave band.
    fc --> current center frequency
    octv --> octave-band (octave-->1, 1/3 octave-->1/3)
    """
    upper = 2 ** ( octv / 2) * fc
    lower = 2 ** (-octv / 2) * fc
    return upper, lower

def OctaveCenterFreqsGen(dx=3, n=39):
    """Produce general center frequencies for octave-band spectra
    dx --> frequency interval spacing (3 for octave, 1 for 1/3 octave)
    n --> number of center freqs to product (starting at dx)
    """
    fc30 = 1000 #Preferred center freq for m=30 is 1000Hz
    m = np.arange(1, n+1) * dx #for n center freqs, multiply 1-->n by dx
    freqs = fc30 * 2 ** (-10 + m/3) #Formula for center freqs

#def OctaveCenterFreqs(narrow, octv=1):
#    """Calculate center frequencies (fc) for octave or 1/3 octave bands.
#    Provide original narrow-band frequency vector to bound octave-band.
#    narrow --> original narrow-band frequencies (provides bounds for octave)
#    octv --> frequency interval spacing (1 for octave, 1/3 for 1/3 octave)
#    """
#    fc30 = 1000 #Preferred center freq for m=30 is 1000Hz
#    freqs = []
#    for i in range(len(narrow)):
#        #current index
#        m = (3 * octv) * (i + 1) #octave, every 3rd, 1/3 octave, every 1
#        freq = fc30 * 2 ** (-10 + m/3) #Formula for center freq
#        if freq > max(narrow):
#            break #quit if current fc is greater than original range
#        if freq >= min(narrow):
#            freqs.append(freq) #if current fc is in original range, save
#    return freqs

def OctaveCenterFreqs(narrow, octv=1):
    """Calculate center frequencies (fc) for octave or 1/3 octave bands.
    Provide original narrow-band frequency vector to bound octave-band.
    narrow --> original narrow-band frequencies (provides bounds for octave)
    octv --> frequency interval spacing (1 for octave, 1/3 for 1/3 octave)
    """
    fc30 = 1000 #Preferred center freq for m=30 is 1000Hz
    freqs = []
    for i in range(len(narrow)):
        #current index
        m = (3 * octv) * (i + 1) #octave, every 3rd, 1/3 octave, every 1
        fc = fc30 * 2 ** (-10 + m/3) #Formula for center freq
        fcu, fcl = OctaveBounds(fc, octv) #upper and lower bounds for fc band
        if fcu > max(narrow):
            break #quit if current fc is greater than original range
        if fcl >= min(narrow):
            freqs.append(fc) #if current fc is in original range, save
    return freqs

def OctaveLp(df, octv=1):
    """Get SPL ( Lp(fc,m) ) for octave-band center frequency.
    Returns octave-band center frequencies and corresponding SPLs
    df --> pandas dataframe containing narrow-band frequencies and SPL
    octv --> octave-band type (octave-->1, 1/3 octave-->1/3)
    """

    #Get Center Frequencies
    fcs = OctaveCenterFreqs(df['freq'], octv)
    Lp_octv = np.zeros(len(fcs))
    for i, fc in enumerate(fcs):
        #Get Upper/Lower center freqency band bounds
        fcu, fcl = OctaveBounds(fc, octv)

        band = df[df['freq'] >= fcl]
        band = band[band['freq'] <= fcu]
        #band = df[df['freq'] >= fcl and df['freq'] <= fcu]

        #SPLs in current octave-band
        Lp = np.array(band['SPL'])
        #Sum 10^(Lp/10) accross current octave-band
        Sum = np.sum( 10 ** (Lp / 10) )
        Lp_octv[i] = 10 * np.log10(Sum)

    return fcs, Lp_octv




def main(source):
    """Perform calculations for frequency data processing
    source --> file name of source sound file
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

    df['SPL'] = SPLi(df['Pa']) #vs time
    powspec['SPL'] = SPLf(Gxx, T) #vs freq

    #print(len(freqs))
    #print(len(df['SPL']))

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
    ### OCTAVE-BAND CONVERSION #########################################
    ####################################################################

    #frq = OctaveCenterFreqs(np.array(powspec['freq']), octv=1)
    #print(frq)

    octv3rd = pd.DataFrame()
    octv3rd['freq'], octv3rd['SPL'] = OctaveLp(powspec, octv=1/3)

    octv = pd.DataFrame()
    octv['freq'], octv['SPL'] = OctaveLp(powspec, octv=1)









    ####################################################################
    ### SAVE DATA ######################################################
    ####################################################################

    #SAVE WAVE SIGNAL DATA
    df = df[['time', 'Pa', 'SPL', 'V']] #reorder
    df.to_csv( '{}/timespec.dat'.format(datadir), sep=' ', index=False ) #save

    #SAVE POWER SPECTRUM DATA
    powspec.to_csv( '{}/freqspec.dat'.format(datadir), sep=' ', index=False )

    #SAVE OCTAVE-BAND DATA
    octv3rd.to_csv( '{}/octv3rd.dat'.format(datadir), sep=' ', index=False)
    octv.to_csv( '{}/octv.dat'.format(datadir), sep=' ', index=False)

    #SAVE SINGLE PARAMETERS
    params = pd.DataFrame()
    params = params.append(pd.Series({'fs' : fs, 'tNwave' : dt_Nwave, }), ignore_index=True)
    params.to_csv( '{}/params.dat'.format(datadir), sep=' ', index=False)


if __name__ == "__main__":


    Source = 'Boom_F1B2_6.wav'

    main(Source)
