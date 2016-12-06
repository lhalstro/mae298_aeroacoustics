"""DATA PROCESSING PROGRAM
Logan Halstrom
MAE 298 AEROACOUSTICS
HOMEWORK 4 - TURBOMACHINERY NOISE
CREATED: 17 NOV 2016
MODIFIY: 06 DEC 2016

DESCRIPTION:  Solve for upstream propagation modes in an axial jet turbine
engine.  Perform sound pressure level analysis for blade self noise caused
by a gust entering the inlet and encountering a fan blade.

NOTE:  Calculations are for a 1/6 scale model of the jet engine and z=0 plane
pressure data is derived from a CFD model
"""

#IMPORT GLOBAL VARIABLES
from hw4_98_globalVars import *

import numpy as np
import pandas as pd

from scipy.special import jn, yn #bessel functions
from scipy.special import jv, yv #bessel functions
from scipy.optimize import fsolve #linear solver
from scipy.optimize import broyden1 #non-linear solver

#CUSTOM PLOTTING PACKAGE
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/Logan/lib/python')
from lutil import df2tex
from lplot import *
from seaborn import color_palette
import seaborn as sns
UseSeaborn('xkcd') #use seaborn plotting features with custom colors
colors = sns.color_palette() #color cycle
markers = bigmarkers         #marker cycle


def AxialEigenvalFunc(x, m, ri, ro):
    """Determinant of axial flow boundary condition eigenfunctions, which is
    equal to zero.  Use scipy solver to get eigenvalues
    x  --> dependent variable (eigenvalue mu in axial flow eigenfunction)
    m  --> Order of bessel function
    ri --> inner radius of jet engine
    ro --> ourter radius of jet engine
    """

    X = x * ro
    Y = x * ri
    return ( 0.5 * (jv(m-1, X) - jv(m+1, X)) * 0.5 * (yv(m-1, Y) - yv(m+1, Y))
           - 0.5 * (jv(m-1, Y) - jv(m+1, Y)) * 0.5 * (yv(m-1, X) - yv(m+1, X))
            )

def AxialEigenfunction(r, ri, m, mu):
    """Axial flow eigenfunction
    r  --> radial location to evaluate eigenfunction at
    ri --> inner radius of axial jet engine
    m  --> circumfrential acoustic mode
    mu --> eigenvalue (dependent on radial mode n)
    """
    X = mu * ri
    return ( jv(m, mu * r) - (0.5 * (jv(m-1, X) - jv(m+1, X)))
                / (0.5 * (yv(m-1, X) - yv(m+1, X))) * yv(m, mu * r) )


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

def CutOnCondition(Kz):
    """Determine if mode is cuton based on wave number
    Must be real to propagate (pos. imag.: decay exponential, neg.: expand)
    Must be negative for upstream propagation
    Kz --> axial wavenumber for current mode
    """
    cuton = 'No'
    #Test if no imaginary part (propagation)
    if Kz.imag == 0:
        #Test if upstream propagation
        if Kz.real < 0:
            cuton = 'Yes'
    return cuton


def GetGammaAmn(m, n, mu, p, r):
    """Get Gamma_mn and Amn from eigenfunction and acoustic pressure
    m     --> current circumfrential mode
    n     --> current radial mode
    mu    --> eigenvalue for current m,n
    p   --> acoustic pressure at z=0 plane
    r   --> radius vector associated with p (goes from Ri --> Ro)
    """
    ri, ro = min(r), max(r) #inner/outer radii of engine

    #Calculate Gamma_mn for non m=n=0 case:
    Gam = (0.5 * (ro ** 2 - m ** 2 / mu ** 2)
            * AxialEigenfunction(ro, ri, m, mu) ** 2
         - 0.5 * (ri ** 2 - m ** 2 / mu ** 2)
            * AxialEigenfunction(ri, ri, m, mu) ** 2 )
    #Calculate Amn
    # Psi = lambda r: AxialEigenfunction(r, ri, m, mu)
    Psi = AxialEigenfunction(r, ri, m, mu)
    Amn = 1 / Gam * np.trapz( p * Psi * r, r)


    return Gam, Amn

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

    # #tweak for round off error
        #     #values for calculating Wmn are very near zero.  small round-off
        #     #errors can cause Wmn=0 and PWL=infinity.
        #     #Change number of sigfigs to get correct Wmn for all cases
    rounds = [] #round eigenvalue to deal with floating point error
    for m in ms:
        rounds.append([None, None, None, None, None])
    rounds[0] = [8, None, None, None, None]

    for i, m in enumerate(ms):

        #guesses for root solutions
        x0s = np.linspace(0.1, 3, 300)
        mus = []
        for j, x0 in enumerate(x0s):

            #try each solution guess to see if it returns a good result
            try:
                #find solution corresponding to current guess
                mu = broyden1( lambda x: AxialEigenvalFunc(x, m, Ri, Ro), x0)
                mu = float(mu)
                #add to solution list if no error message
                mus.append( mu )
            except Exception:
                #Skip any solution errors
                pass

        #GET ONLY UNIQUE VALUES OF MU
        df = pd.DataFrame({'long' : list(mus), 'short' : list(mus)})
        df = df.round({'short': 6}) #shorten values to find unique ones
        df = df.drop_duplicates(subset='short') #drop duplicates in short vals
        df = df.sort_values('long') #sort eigen values from least to greatest
        df = df.reset_index() #reset indicies to same as n

        #Asign desired number of eigen values to solution dataframe
        mus = df['long'][:Nn]
        #Deal with floating point error
        for j, mu in enumerate(mus):
            if rounds[i][j] != None:
                mus[j] =  round(mu, rounds[i][j])
        eigenvals[m] = mus
        # eigenvals = eigenvals.round({m : 8})

    #SAVE EIGENVALUES
    #columns are m, rows are n
    eigenvals.to_csv( '{}/eigenvalues.dat'.format(datadir), sep=' ',
                        index=True )
    df2tex(eigenvals, '{}/eigenvalues'.format(datadir), dec=sigfigs)



    ####################################################################
    ### PROB 2 - PLOT EIGENFUNCTIONS ###################################
    ####################################################################

    R = np.linspace(Ri, Ro, 201) #Radial vector in engine
    line = ['-', '--', '-.', ':']

    fig, ax = plt.subplots(len(ms), sharex=True, figsize=[3, 12])
    for i, m in enumerate(ms):
        for j, n in enumerate(ns):
            ax[i].plot(R, AxialEigenfunction(R, Ri, m, eigenvals[m][j]),
                        label='n={}'.format(j)
                        )
            ax[i].set_ylabel('$\\Psi_{{{},n}}$'.format(m))
            ax[i].set_xlim([Ri, Ro])
    ax[i].set_xlabel('Radial Location [in]') #label x-axis on last subplot

    savename = '{}/2_eigenfunctions_subplot.{}'.format(picdir, pictype)
    SavePlot(savename)



    ###############

    _,ax = PlotStart(None, 'Radial Location [in]',
                            'Eigenfunction ($\\Psi_{{m,n}}$)', figsize=[6, 6])

    ax.plot([Ri, Ro], [0, 0], color='black', linestyle='--', linewidth=2, alpha=0.7) #zero
    mhandles = []
    mlabels = []
    nhandles = []
    nlabels = []
    labels = []
    for i, m in enumerate(ms):
        for j, n in enumerate(ns):
            h, = ax.plot(R, AxialEigenfunction(R, Ri, m, eigenvals[m][j]),
                        label='n={}'.format(j), color=colors[j],
                        linestyle=line[i],
                        )
            if j == 0:
                mhandles.append(h)
                mlabels.append(str(m))
            if i == 0:
                nhandles.append(h)
                nlabels.append(str(n))
            labels.append(j)
    ax.set_xlim([Ri, Ro])

    leg1 = PlotLegendLabels(ax, mhandles, mlabels, loc='upper left', title='m')
    leg2 = PlotLegendLabels(ax, nhandles, nlabels, loc='lower left', title='n')
    plt.gca().add_artist(leg1)

    savename = '{}/2_eigenfunctions.{}'.format(picdir, pictype)
    SavePlot(savename)


    ####################################################################
    ### PROB 3 - WAVE NUMBER/CUT-OFF CONDITION #########################
    ####################################################################


    wavenums = eigenvals.copy()
    wavenums = wavenums.applymap(lambda x: AxialWavenumber(x, omega, a, M))

    #SAVE EIGENVALUES
    #columns are m, rows are n
    wavenums.to_csv( '{}/wavenumbers.dat'.format(datadir), sep=' ',
                        index=True )

    #DETERMINE CUTON CONDITION
    curNs = [0, 1, 2]
    #sort by m
    for m in ms:
        Kzs = []
        for n in curNs:
            Kzs.append(wavenums[m][n])

        df = pd.DataFrame({'m' : [m for x in curNs], 'n' : curNs, 'Kz' : Kzs})
        df['Cut-on'] = [CutOnCondition(Kz) for Kz in df['Kz']]
        #save wavenumber and cuton table
        df = df[['m', 'n', 'Kz', 'Cut-on']]
        df2tex(df, '{}/wavenumbers_m{}'.format(datadir, m), dec=sigfigs)


    ####################################################################
    ### PROB 4 - SPL OF REAL PRESSURE FIELD ############################
    ####################################################################

    #READ EXPERIMENTAL PRESSURE DISTRIBUTION AT Z=0 PLANE
        #Columns: Radius [in], Real pressure [psi], Imaginary pressure [psi]
    df = pd.read_csv('{}/pressure_input.dat'.format(datadir), sep='\t',
                    names=['R', 'pRe', 'pIm'])
    #Create complex values of acoustic pressure
    df['p'] = df['pRe'] + df['pIm'] * 1j

    m, n = 18, 0
    PWLs = []

    for i, n in enumerate([0, 1, 2]):
        mu = eigenvals[m][n] #eigenvalue for curent (m,n)
        Gam, Amn = GetGammaAmn(m, n, mu, df['p'], df['R'])

        #MODAL POWER
        Kz = wavenums[m][n]
        frac = Kz / (omega / a - Kz * M)
        Wmn = np.pi / (rho * a) * Gam * Amn * np.conj(Amn) * (
                (1 + M ** 2) * frac.real + M * (abs(frac) ** 2 + 1) )

        #SOUND POWER LEVEL
        PWL = 10 * np.log10( abs(Wmn) ) - 10* np.log10(7.3756e-13)

        PWLs.append(PWL)

    #SAVE POWERLEVELS
    PWLs = pd.DataFrame({m : PWLs})
    PWLs.to_csv( '{}/powerlevels.dat'.format(datadir), sep=' ', index=True )
    df2tex(PWLs, '{}/powerlevels'.format(datadir), dec=sigfigs)


if __name__ == "__main__":


    main()
