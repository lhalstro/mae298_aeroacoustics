"""WRAPPER PROGRAM
Logan Halstrom
MAE 298 AEROACOUSTICS
HOMEWORK 4 - TURBOMACHINERY NOISE
CREATED: 17 NOV 2016
MODIFIY: 17 NOV 2016

DESCRIPTION: Run data processing and plotting programs in order.
Enter common inputs in this script.
"""

#IMPORT GLOBAL VARIABLES AND PROGRAMS TO WRAP
from hw1_98_globalVars import *
import hw1_00_process as process
import hw1_01_plot as plot


def main(source):
    """Wrapper script for project.  First compute and save relevant data.
    Then load and plot data.
    source --> path to file containing source pressure data at z=0 plane
    """

    print('\nProcessing Data')
    process.main(source)
    print('\nPlotting Data')
    plot.main()


if __name__ == "__main__":

    Source = '{}/pressure_input.dat'.format(datadir)

    main(Source, InletSpacing)
