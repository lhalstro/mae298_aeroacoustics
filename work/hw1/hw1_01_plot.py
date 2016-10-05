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
    df = pd.read_csv('{}/signal.dat'.format(datadir), sep=' ' )
    powspec = pd.read_csv('{}/powspec.dat'.format(datadir), sep=' ' )


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

    savename = '{}/1_3_SPL.{}'.format(picdir, pictype)
    SavePlot(savename)





if __name__ == "__main__":




    main()
