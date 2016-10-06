"""HW1 - DATA PLOTTING
Logan Halstrom
MAE 298 AEROACOUSTICS
HOMEWORK 1 - SIGNAL PROCESSING
CREATED: 04 OCT 2016
MODIFIY: 04 OCT 2016

DESCRIPTION: Plot processed signal of sonic boom.
narrow-band spectrum
single-side spectral density
SPL
octave bands
overall SPL
"""

#IMPORT GLOBAL VARIABLES
from hw1_98_globalVars import *

import numpy as np
import pandas as pd

import os

#CUSTOM PLOTTING PACKAGE
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/Logan/lib/python')
from lplot import *
from seaborn import color_palette
import seaborn as sns
UseSeaborn('xkcd') #use seaborn plotting features with custom colors
colors = sns.color_palette() #color cycle
markers = bigmarkers         #marker cycle

MarkerWidth = 2.25

def main():
    """input description
    """

    #Make Plot Output Directory
    MakeOutputDir(picdir)

    #LOAD DATA
    df = pd.read_csv('{}/timespec.dat'.format(datadir), sep=' ' )
    powspec = pd.read_csv('{}/freqspec.dat'.format(datadir), sep=' ' )
    params = pd.read_csv('{}/params.dat'.format(datadir), sep=' ' )
    octv3rd = pd.read_csv('{}/octv3rd.dat'.format(datadir), sep=' ' )
    octv = pd.read_csv('{}/octv.dat'.format(datadir), sep=' ' )


    ####################################################################
    ### 1.1 PLOT PRESSURE WAVE #########################################
    ####################################################################

    _,ax = PlotStart(None, 'Time (s)', 'Pressure (MPa)', figsize=[6, 6])
    #Hollow Marker Plot
    ax.plot(df['time'], df['Pa'] / 10e6,
            #label=lbl, color=clr,
            # linewidth=0,
            marker=markers[0], markevery=500,
            markeredgecolor=colors[0], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )

    savename = '{}/1_1_Pressure.{}'.format(picdir, pictype)
    SavePlot(savename)

    ####################################################################
    ### 1.2 PLOT POWER SPECTRUM DENSITY ################################
    ####################################################################

    _,ax = PlotStart(None, 'Frequency (Hz)', 'Power ($G_{xx}$)', figsize=[6, 6])
    #Hollow Marker Plot
    ax.plot(powspec['freq'], powspec['Gxx'],
            #label=lbl, color=clr,
            # linewidth=0,
            marker=markers[0], markevery=1,
            markeredgecolor=colors[0], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )
    plt.xlim([0,50])

    savename = '{}/1_2_PowerSpec.{}'.format(picdir, pictype)
    SavePlot(savename)

    ####################################################################
    ### 1.3 PLOT SOUND PRESSURE LEVEL ##################################
    ####################################################################

    _,ax = PlotStart(None, 'Time (s)', 'SPL (dB)', figsize=[6, 6])
    #Hollow Marker Plot
    ax.plot(df['time'], df['SPL'],
            #label=lbl, color=clr,
            # linewidth=0,
            marker=markers[0], markevery=500,
            markeredgecolor=colors[0], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )

    savename = '{}/1_3_SPLt.{}'.format(picdir, pictype)
    SavePlot(savename)

    _,ax = PlotStart(None, 'Frequency (Hz)', 'SPL (dB)', figsize=[6, 6])
    #Hollow Marker Plot
    ax.plot(powspec['freq'], powspec['SPL'],
            #label=lbl, color=clr,
            # linewidth=0,
            marker=markers[0], markevery=500,
            markeredgecolor=colors[0], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )
    ax.set_xscale('log')

    savename = '{}/1_3_SPLf.{}'.format(picdir, pictype)
    SavePlot(savename)

    ####################################################################
    ### 2.1 PLOT 1/3 OCTAVE-BAND SPL ###################################
    ####################################################################

    _,ax = PlotStart(None, 'Frequency (Hz)', 'SPL (dB)', figsize=[6, 6])
    #Hollow Marker Plot
    ax.plot(octv3rd['freq'], octv3rd['SPL'],
            #label=lbl, color=clr,
            # linewidth=0,
            marker=markers[0], markevery=500,
            markeredgecolor=colors[0], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )
    ax.set_xscale('log')

    savename = '{}/2_1_SPLf_octv3rd.{}'.format(picdir, pictype)
    SavePlot(savename)

    ####################################################################
    ### 2.2 PLOT OCTAVE-BAND SPL #######################################
    ####################################################################

    _,ax = PlotStart(None, 'Frequency (Hz)', 'SPL (dB)', figsize=[6, 6])
    #Hollow Marker Plot
    ax.plot(octv['freq'], octv['SPL'],
            #label=lbl, color=clr,
            # linewidth=0,
            marker=markers[0], markevery=500,
            markeredgecolor=colors[0], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )
    ax.set_xscale('log')

    savename = '{}/2_2_SPLf_octv.{}'.format(picdir, pictype)
    SavePlot(savename)

    ####################################################################
    ### PLOT ALL OCTAVE-BAND SPL #######################################
    ####################################################################

    _,ax = PlotStart(None, 'Frequency (Hz)', 'SPL (dB)', figsize=[6, 6])

    #dats = [powspec, octv3rd, octv]
    #labels = ['narrow', 'octave', '$\\frac{1}{3}$octave']
    #mkevery = [25, 1, 1]

    #for dat, lbl, clr, mkr, mkev in zip(dats, labels, colors, markers, mkevery):
    #    ax.plot(dat['freq'], dat['SPL'], label=lbl,
    #            #color=clr,
    #            # linewidth=0,
    #            marker=mkr, markevery=mkev,
    #            markeredgecolor=clr, markeredgewidth=MarkerWidth,
    #            markerfacecolor="None",
    #            )

    #Plot Overall SPL as Horizontal Line
    xmin = min(powspec['freq'])
    xmax = max(powspec['freq'])
    ax.plot([xmin, xmax], [params['SPL_overall'], params['SPL_overall']],
            label='Overall', color='black', linestyle='--')

    i = 1
    #Plot Octave-Band
    ax.plot(octv['freq'], octv['SPL'], label='oct',
            color=colors[i],
            marker=markers[i], markevery=1,
            markeredgecolor=colors[i], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )

    i += 1
    #Plot 1/3 Octave-Band
    ax.plot(octv3rd['freq'], octv3rd['SPL'], label='$\\frac{1}{3}$oct',
            color=colors[i],
            marker=markers[i], markevery=1,
            markeredgecolor=colors[i], markeredgewidth=MarkerWidth,
            markerfacecolor="None",
            )

    #i += 1
    i = 0
    #Plot Narrow-Band
    ax.plot(powspec['freq'], powspec['SPL'], label='nar',
            color=colors[i],
            )

    ax.set_xscale('log')
    ax.set_xlim([xmin, xmax])
    PlotLegend(ax)

    savename = '{}/2_SPLf_all.{}'.format(picdir, pictype)
    SavePlot(savename)



if __name__ == "__main__":




    main()
