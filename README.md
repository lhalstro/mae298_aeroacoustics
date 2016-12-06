# Aeroacoustics

##### MAE 298 - Fall Quarter 2016<br>UC Davis<br>Prof: Dr. Seongkyu Lee

A general introduciton to aeroacoustics beginning with wave propagation and perception of noise and then advancing to analytic techniques including Lighthill's acoustic analogy and Lilly's equation.  The course is concluded with the Ffwocks-Williams Hawkings equations and axial flow relations.

### Repo Organization
* ClassNotes
    * Scans of lecture notes for the course
    * Contains first day presentation, which acts as course syllabus
* hwX
    * Directories for each homework
    * Contains code and writeup
    * Contains scaned corrected versions of student submissions
    * See below for details
* Midterm
    * Solution of miterm exam
    * Scan of student submission with corrections
* FinalProject
    * Literature study on the topic of Prediction Methods in Aeroacoustics of Launch Vehicles
    * 'references' directory is not tracked by GIT

### Student Work
1. Acoustic Signal Processing and Sound Pressure Level/Octave Band Analysis
  * Analyzes recorded sonic boom pressure signal
  * Converts pressure signal to FFT frequency spectrum
  * Calculate Sound Pressure Level of signal in 1/3 Octave and Octave bands
2. Lilley's Equation Solution and Application to Jet Noise
  * Derivation of Lilley's equation for jet noise
  * Application to flow regions inside and outside jet potential core
  * Discussion of necessary boundary conditions for matching two flows
3. Generalized Differentiation and Farassat's Formulation
  * Generalized differentiation of the wave equation to derive Kirchoff formulat for moving boties
  * Derivation of Farassat Formulation 1A for Loading Noise
4. Turbomachinery Noise
  * Solve for upstream propagation modes in an axial jet turbine engine
  * Perform sound pressure level analysis for blade self noise caused by a gust entering the inlet and encountering a fan blade
5. Final Project Literature Study
  * Review of Aeroacoustic Prediction Techniques for Launch Vehicles
  * Experimental Testing
    * free-field, reverberant chamber, progressive wave testing
    * Flight testing, microphone phased-array
  * Wind Tunnel Testing
    * Sub-scale, high-fidelity geometry model testing at simulated flight conditions
    * Heated helium abort motor plume simulation
  * Computational Methods
    * Computational Fluid Dynamics
        * Accurate unsteady surface static pressure
        * Inaccurate aeroacoustic pressure fluctuations
    * Computational AeroAcoustics
        * CFD/CAA coupling
        * Needs improvement
