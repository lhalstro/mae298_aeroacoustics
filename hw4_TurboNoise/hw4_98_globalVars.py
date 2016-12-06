"""GLOBAL VARIABLES STORAGE
Logan Halstrom
MAE 298 AEROACOUSTICS
HOMEWORK 4 - TURBOMACHINERY NOISE
CREATED: 17 NOV 2016
MODIFIY: 06 DEC 2016

DESCRIPTION: Provide global variables for all scripts including wrapper.
"""

import numpy as np

#INPUTS

#CFD Engine Model Parameters
scale = 1/6 #CFD engine model scale
Nblade = 18 #Number of blades in engine fan
Nvane  = 1  #Uneven gust loading modeled as single vane
Ro = 13     #Model engine (outer) radius [in]
Ri = 3      #Model hub (inner) radius [in]
Linlet = 4  #Distace between inlet and z=0 plane [in]

#Flow Parameters
M = 0.525        #Engine flow mach number
a = 13503.937009 #Speed of Sound [in/s]
rho = 1.4988e-5 #density [slug/in^3]
RPM = 8326.3042  #Fan RPM
omega = RPM * (2 * np.pi / 60) * Nblade #angular frequency [rad/s]

#DATA OVERWRITE SWITCHES
datoverwrite = 1 #Overwrite data = 1

#LOAD/SAVE DIRECTORIES
datadir = 'Data'    #Source and processed data storage directory
savedir = 'Results' #

picdir = 'Plots' #Plot storage directory
pictype = 'png'      #Plot save filetype
# pictype = 'pdf'      #Plot save filetype

sigfigs = 4 #number of sig figs to save in data files
