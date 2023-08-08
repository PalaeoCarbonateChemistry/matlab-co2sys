ChangeLog


1) Can accept input parameters of [CO3], [HCO3], and [CO2], and propagate their uncertainties
2) Includes NH3 and HS as alkalinity contributors, and propagates their uncertainties
3) Uses separate inputs to specify choices for characterizations of K1K2, KSO4, KF, and TB
4) Does not evaluate input parameters equal to -999 or NaN
5) Exits pH iteration loops that do not converge and indicates where a problem occurred
6) Provides exactly identical pH results for a given input line, no matter the other lines of input parameters (this is not necessarily the case for prior versions of CO2SYS.m)
7) Uses an updated definition of the ideal gas constant (https://physics.nist.gov/cgi-bin/cuu/Value?r)
8) Fixes bugs in CO2SYS.m Revelle factor calculation and derivnum.m output conditions
9) Includes K1 and K2 constants defined by Sulpis et al. (2020), K2 constant defined by Schockman and Byrne (2021), KF constant defined by Perez and Fraga (1987), and KSO4 constant of Waters and Millero (2013)
10) Determines initial pH in iterative solvers using the approach of Munhoven (2013), detailed further in Humphreys et al. (submitted), rather than simply using an initial estimate of 8.0 each time.
11) Obtains free scale pH properly within iterative pH solvers no matter the input scale, rather than making the simplification that input pH is always on the total scale.
12) Includes substrate-inhibitor ratio (Bach, 2015) as an output argument from CO2SYS.
13) Calculates uncertainties at output conditions that are associated with equilibrium constants with respect to equilibrium constants at output conditions, rather than input conditions as previously. This essentially assume pK uncertainty is constant regardless of temperature and pressure.
14) Calculates derivatives and errors for the Revelle factor.
15) errors.m includes optional calcium concentration uncertainty input as discussed in Dillon et al. (2020)
16) Option added for pressure corrections to K0 and fugacity factor.


version 2.0.5, 2018-11-23
-----------------------------------
	- fixed error in eBt x deriv (array) - Thanks to Anna Conchon for bug report

version 2.0.4, 2018-11-10
-----------------------------------
	- updated README; Orr et al. (2018) published
	- changed defaults of epK and eBt, consistent w/ published version of Orr et al. (2018)

version 2.0.3, 2018-06-04
-----------------------------------
	- changed default uncertainties for epK and eBt for consistency with Orr et al. (2018)
	- new Jupyter notebooks added to demo use of errors, derivnum, CO2SYS routines
	- corrected derivnum.m units of PAR1ref & PAR2ref [H+]
	- fixed eBt error calc in errors.m

version 2.0.2, 2017-10-17
-----------------------------------
	- corrections so that errors.m works when there a vectors in input
	- changed Octave enhancemnts (!, size, resize) to be MATLAB compatible (~, ...) 

version 2.0.1, 2017-10-11
-----------------------------------
	- new routines to support use of TEOS-10 standards (Conserrative temperature, Aboslute salinity)
	- added uncertainty in total boron (eBt) as a separate argument in errors.m
	- new and improved examples

version 2.0, 2016-08-15
-----------------------------------
	- errors.m: new routine added to propagate uncertainties
	- derivnum.m: new routine added to compute partial derivatives, needed for error propagation
	- CO2SYS.m: slightly modified routine 
		* slight changes to allow error propagation
		* new option to choose K1 & K2 from Waters et al. (2014): fixes inconsistencies with Millero (2010) identified by Orr et al. (2015)

version 1.1, 2016-12-20
-----------------------------------
	- Original version of CO2SYS.m from van Heuven et al. (2011) 
	- downloaded on 20 Dec 2016 from CDIAC: http://cdiac.ornl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/


REFERENCES

Millero, F. J.: Carbonate constants for estuarine waters, Mar. Freshwater Res., 61, 139–142, doi:10.1071/MF09254, 2010.

Orr, J. C. and Epitalon, J.-M.: Improved routines to model the ocean carbonate system: mocsy 2.0, Geosci. Model Dev., 8, 485–499,
doi:10.5194/gmd–8–485–2015, 2015.

Orr, J. C., Gattuso, J.-P., and Epitalon, J.-M.: Comparison of ten packages that compute ocean carbonate chemistry, Biogeosciences, 12,
35 1483–1510, doi:10.5194/bg–12–1483–2015, 2015.

van Heuven, S., D. Pierrot, J.W.B. Rae, E. Lewis, and D.W.R. Wallace. 2011. MATLAB Program Developed for CO2 System
Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, 
U.S. Department of Energy, Oak Ridge, Tennessee. doi: 10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1

Waters, J., Millero, F., and Woosley, R.: Corrigendum to “The free proton concentration scale for seawater pH”,[MARCHE: 149 (2013)
10 8-22], Mar. Chem., 165, 66–67, 2014.
