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
from scipy.optimize import fsolve #linear solver
from scipy.optimize import broyden1 #non-linear solver



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


def AxialEigenvalFunc(x, m, ri, ro):
    """Determinant of axial flow boundary condition eigenfunctions, which is
    equal to zero.  Use scipy solver to get eigenvalues
    x  --> dependent variable (eigenvalue mu in axial flow eigenfunction)
    m  --> Order of bessel function
    ri --> inner radius of jet engine
    ro --> ourter radius of jet engine
    """
    # return (0.5 * (jn(m-1, x * ro) - jn(m+1, x * ro))
    #                   * 0.5 * (yn(m-1, x * ri) - yn(m+1, x * ri))
    #                   - 0.5 * (jn(m-1, x * ri) - jn(m+1, x * ri))
    #                   * 0.5 * (yn(m-1, x * ro) - yn(m+1, x * ro))
    #         )

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

def AxialEigenfunctionLambda(ri, m, mu):
    """Axial flow eigenfunction, with lambda function radial input for
    integration
    ri --> inner radius of axial jet engine
    m  --> circumfrential acoustic mode
    mu --> eigenvalue (dependent on radial mode n)
    """
    return lambda r: jn(m, mu * r) - jn(m, mu * ri) / yn(m, mu * ri) * yn(m, mu * r)

def AxialWavenumber(mu, omega, c, M):
    """Calculate z-direction wavenumber (Kz+) for sound mode in axial flow
    mu --> eigenvalue of current mode
    omega --> angular frequency of flow
    c --> speed of sound
    M --> mach number of flow
    """
    K = omega / c #flow wave number
    Kz = (-M + np.emath.sqrt(1 - (1 - M ** 2) * (mu / K) ** 2 )) / (1 - M ** 2) * K
    return Kz



def main():
    """Perform calculations for frequency data processing
    """

    ####################################################################
    ### PROB 1 - FIND EIGENVALUES ######################################
    ####################################################################


    eigenvals = pd.DataFrame() #solutions for eigenvalues

    ms = [18, 17, 16, 15] #circumfrential modes to solve for
    Nn = 5 #number of radial modes to solve for
    ns = range(Nn) #find first five radial modes

    for m in ms:

        #guesses for root solutions
        x0s = np.linspace(0.1,3, 30)
        mus = []
        for x0 in x0s:

            #try each solution guess to see if it returns a good result
            try:
                #find solution corresponding to current guess
                mu = broyden1( lambda x: AxialEigenvalFunc(x, m, Ri, Ro), x0)
                #add to solution list if no error message
                    #convert to float and round of floating point error
                mus.append( round(float(mu), 6) )
            except Exception:
                #Skip any solution errors
                pass

        #Get only unique solutions with 'set'
            #convert to numpy array and sort least to greatest
        mus = np.sort(np.array( list(set(mus)) ), axis=None)

        eigenvals[m] = mus[:Nn]

    #SAVE EIGENVALUES
    #columns are m, rows are n
    eigenvals.to_csv( '{}/eigenvalues.dat'.format(datadir), sep=' ',
                        index=True )



    ####################################################################
    ### PROB 2 - PLOT EIGENFUNCTIONS ###################################
    ####################################################################


    # R = np.linspace(Ri, Ro, 101)

    # plt.figure()
    # plt.plot(R, AxialEigenfunction(R, Ri, m, mu) )

    # plt.xlim([Ri, Ro])
    # plt.show()




    R = np.linspace(Ri, Ro, 101) #Radial vector in engine

    fig, ax = plt.subplots(len(ms), sharex=True, figsize=[5, 10])
    for i, m in enumerate(ms):
        for j, n in enumerate(ns):
            ax[i].plot(R, AxialEigenfunction(R, Ri, m, eigenvals[m][j]),
                        label='n={}'.format(j)
                        )
            ax[i].set_ylabel('$\\Psi_{{{},n}}$'.format(m))
            ax[i].set_xlim([Ri, Ro])
    ax[i].set_xlabel('Radial Location [in]') #label x-axis on last subplot

    plt.show()



    ###############

    _,ax = PlotStart(None, 'Radial Location [in]', '$\\Psi_{{m,n}}$', figsize=[6, 6])

    for i, m in enumerate(ms):

        for j, n in enumerate(ns):

            ax.plot(R, AxialEigenfunction(R, Ri, m, eigenvals[m][j]),
                        label='n={}'.format(j), color=colors[j],
                        )
    # ax.set_ylabel('$\\Psi_{{{}n}}$'.format(m))
    # ax.set_xlabel('Radius [in]') #label x-axis on last subplot
    ax.set_xlim([Ri, Ro])

    plt.show()


    ####################################################################
    ### PROB 3 - WAVE NUMBER/CUT-OFF CONDITION #########################
    ####################################################################


    wavenums = eigenvals.copy()

    # wavenums = AxialWavenumber(eigenvals, omega, a, M)
    wavenums = wavenums.applymap(lambda x: AxialWavenumber(x, omega, a, M))

    #SAVE EIGENVALUES
    #columns are m, rows are n
    wavenums.to_csv( '{}/wavenumbers.dat'.format(datadir), sep=' ',
                        index=True )

    #Cut-on if real and negative




    ####################################################################
    ### PROB 4 - SPL OF REAL PRESSURE FIELD ############################
    ####################################################################

    #READ EXPERIMENTAL PRESSURE DISTRIBUTION AT Z=0 PLANE
        #Columns: Radius [in], Real pressure [psi], Imaginary pressure [psi]
    df = pd.read_csv('{}/pressure_input.dat'.format(datadir), sep='\t',
                    names=['R', 'pRe', 'pIm'])

    print(df)










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


    main()
