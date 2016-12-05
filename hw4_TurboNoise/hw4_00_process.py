"""DATA PROCESSING PROGRAM
Logan Halstrom
MAE 298 AEROACOUSTICS
HOMEWORK 4 - TURBOMACHINERY NOISE
CREATED: 17 NOV 2016
MODIFIY: 04 DEC 2016

DESCRIPTION:

NOTE:
"""

#IMPORT GLOBAL VARIABLES
from hw4_98_globalVars import *

import numpy as np
import pandas as pd

from scipy.special import jn, yn #bessel functions
from scipy.optimize import fsolve




def AxialEigenvalFunc(x, *args):
    """Determinant of axial flow boundary condition eigenfunctions, which is
    equal to zero.  Use scipy solver to get eigenvalues
    x    --> dependent variable (eigenvalue mu in axial flow eigenfunction)
    args --> tuple of other arguments: (m, n, ri, ro)
        m --> Order of bessel function
        ri --> inner radius of jet engine
        ro --> ourter radius of jet engine
    """
    m, ri, ro = args

    X = x * ro
    Y = x * ri

    return ( 0.5 * (jn(m-1, X) - jn(m+1, X)) * 0.5 * (yn(m-1, Y) - yn(m+1, Y))
           - 0.5 * (jn(m-1, Y) - jn(m+1, Y)) * 0.5 * (yn(m-1, X) - yn(m+1, X))
            )

def AxialEigenfunction(r, ri, m, mu):
    """Axial flow eigenfunction
    r  --> radial location to evaluate eigenfunction at
    ri --> inner radius of axial jet engine
    m  --> circumfrential acoustic mode
    mu --> eigenvalue (dependent on radial mode n)
    """
    return jn(m, mu * r) - jn(m, mu * ri) / yn(m, mu * ri) * yn(m, mu * r)



def main(source):
    """Perform calculations for frequency data processing
    source --> file name of source sound file
    """

    ####################################################################
    ### PROB 1 - FIND EIGENVALUES ######################################
    ####################################################################


    m = 15
    x0 = 1
    mu = fsolve(AxialEigenvalFunc, x0, args=(m, Ri, Ro) )
    print(mu)

    # print(AxialEigenvalFunc(1, m, Ri, Ro))


    ####################################################################
    ### PROB 2 - PLOT EIGENFUNCTIONS ###################################
    ####################################################################

    import matplotlib.pyplot as plt

    R = np.linspace(Ri, Ro, 101)

    plt.figure()
    plt.plot(R, AxialEigenfunction(R, Ri, m, mu))
    plt.show()











    # ####################################################################
    # ### READ SOUND FILE ################################################
    # ####################################################################

    # df = pd.DataFrame() #Stores signal data

    # #Read source frequency (fs) and signal in volts normallized between -1&1
    # fs, df['V'] = ReadWavNorm( '{}/{}'.format(datadir, source) ) #Like matlab

    # #Convert to pascals
    # df['Pa'] = df['V'] * volt2pasc



    # print('\nNum. of Points, Narrow-band:'    , len(df))
    # print(  'Num. of Points, 1/3 Octave-band:', len(octv3rd))
    # print(  'Num. of Points, Octave-band:'    , len(octv))

    ####################################################################
    ### SAVE DATA ######################################################
    ####################################################################

    # #SAVE WAVE SIGNAL DATA
    # df = df[['time', 'Pa', 'SPL', 'V']] #reorder
    # df.to_csv( '{}/timespec.dat'.format(datadir), sep=' ', index=False ) #save

    # #SAVE POWER SPECTRUM DATA
    # powspec.to_csv( '{}/freqspec.dat'.format(datadir), sep=' ', index=False )

    # #SAVE OCTAVE-BAND DATA
    # octv3rd.to_csv( '{}/octv3rd.dat'.format(datadir), sep=' ', index=False)
    # octv.to_csv( '{}/octv.dat'.format(datadir), sep=' ', index=False)

    # #SAVE SINGLE PARAMETERS
    # params = pd.DataFrame()
    # params = params.append(pd.Series(
    #     {'fs' : fs, 'SPL_overall' : Lp_overall,
    #      'Pmax' : Pmax, 'tNwave' : dt_Nwave,
    #      'ti' : ti, 'Pi' : Pi, 'tf' : tf, 'Pf' : Pf}
    #     ), ignore_index=True)
    # params.to_csv( '{}/params.dat'.format(datadir), sep=' ', index=False)


if __name__ == "__main__":


    Source = 'Boom_F1B2_6.wav'

    main(Source)
