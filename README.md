[![Run small verification tests](https://github.com/MUADh3i9yL/matlab-co2sys-traditional/actions/workflows/small_tests.yml/badge.svg)](https://github.com/MUADh3i9yL/matlab-co2sys-traditional/actions/workflows/small_tests.yml)

![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fgist.githubusercontent.com%2Frossidae%2Fc68eb447f90f281a543bca7ab1d7a56a%2Fraw%2F9ac172abdacff106379815d86e31720f80dadf57%2Fco2sys-performance-metrics.json&query=%24.iterative&suffix=%20%CE%BCsec%2Fpoint&label=Iterative%20speed&color=%23963023)


![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fgist.githubusercontent.com%2Frossidae%2Fc68eb447f90f281a543bca7ab1d7a56a%2Fraw%2F9ac172abdacff106379815d86e31720f80dadf57%2Fco2sys-performance-metrics.json&query=%24.array&suffix=%20%CE%BCsec%2Fpoint&label=Vectorised%20speed)


CITATION:
van Heuven, S., D. Pierrot, J.W.B. Rae, E. Lewis, and D.W.R. Wallace. 2011.
MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b.
Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S.
Department of Energy, Oak Ridge, Tennessee. doi: 10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1


CO2SYS version 1.1 (Sept 2011)

ABOUT CO2SYS:
This is a MATLAB-version of the original CO2SYS for DOS. CO2SYS calculates and returns a detailed 
state of the carbonate system of oceanographic water samples, if supplied with enough input. Use 
this function as you would use any other Matlab inline funtion, i.e., a=func(b,c). For extended 
details on using the function, please refer to the enclosed help by typing "help CO2SYS" in Matlab. 
For details on the internal workings of the function, please refer to the original  publication of 
Lewis and Wallace at http://cdiac.ess-dive.lbl.gov/oceans/co2rprt.html. Note that this function allows 
for the input of vectors. This means that you can calculate many samples at once. Each of these 
samples can be processed with individual salinities, temperatures, pH scales, dissociation constants, etc.

HISTORY:
Original version for DOS was written by Lewis and Wallace. This was converted to MATLAB by Denis Pierrot 
at CIMAS, University of Miami, Miami, Florida. Vectorization, internal refinements and speed improvements 
were added by Steven van Heuven, University of Groningen, The Netherlands. Although functionality has 
been added, the output of the function has not changed. All versions of co2sys available at CDIAC 
(DOS, Excel for WINDOWS, Excel for MAC, MATLAB) should yield (near-) identical results when supplied 
with identical input. If you discover that they don't or you have a more general bug report, please 
email me, Denis Pierrot or Alex Kozyr (svheuven@gmail.com, Denis.Pierrot@noaa.gov, kozyra@ornl.gov). 

INSTALLING:
Download the m-file "CO2SYS.m" and, optionally, the two examples. Place the file(s) in a location that 
Matlab can see (or add the location of the file(s) to Matlab's search path).
Run either of the examples in Matlab, or start using the main routine straight away.
=======
**CITATION**

- If you use any CO2SYS related software, please cite the original work by Lewis and Wallace (1998).
- If you use CO2SYS.m, please cite van Heuven et al (2011).
- If you use errors.m or derivnum.m, please cite Orr et al. (2018).

**CO2SYS-MATLAB versions**

- 1.1   (Sept 2011): van Heuven et al. (2011) 
- 2.0   (20 Dec 2016): includes uncertainty propagation
- 2.0.1 (11 Oct 2017): supports TEOS-10 standards (conservattive temperature, absolute salinity)
- 2.0.2 (17 Oct 2017): Octave enhancements changed to be MATLAB compatible
- 2.0.3 (4 Jun 2018): examples added as Jupyter notebooks
- 2.0.4 (10 Nov 2018): defaults for standard uncertainties in constants (epK vector and eBt) made consistent with Orr et al. (2018), i.e., final published version 
- 2.0.5 (23 Nov 2018): fixed bug in eBt propagation to deriv array (thanks A. Cochon)
- 2.1   (29 Jun 2020): fixed bug in derivnum affecting OUT results (linked to TEMPOUT); masked derivs of input vars in derivnum

**ABOUT CO2SYS**

Here you will find a MATLAB-version of CO2SYS, originally written for
DOS. CO2SYS calculates and returns a detailed state of the carbonate system for
oceanographic water samples, if supplied with sufficient input.  Use the CO2SYS
function as you would use any other MATLAB inline function, i.e.,
a=func(b,c). For much detail about how to use CO2SYS, simply type "help CO2SYS"
in MATLAB.  The help function also works for the two new uncertainty propagation
routines (errors and derivnum).  For details on the internal workings of CO2SYS,
please refer to the original publication (Lewis and Wallace, 1998) available at
http://cdiac.ornl.gov/oceans/co2rprt.html.  Since CO2SYS and the two new
routines each allow input of vectors, with just one call they can process many
samples.  Each sample may have a different salinity, temperature, pH scale,
dissociation constants, etc.

**HISTORY**

The original version for DOS was written by Lewis and Wallace
(1998). That was translated to MATLAB by Denis Pierrot at CIMAS,
University of Miami, Miami, Florida. Then that code was vectorized,
refined, and optimized for computational speed by Steven van Heuven,
University of Groningen, The Netherlands. Although functionality was
added, the output of the CO2SYS function has not changed in form. All
versions of CO2SYS that are available at CDIAC (DOS, Excel, MATLAB)
should produce nearly identical results when supplied with identical
input. Indeed, close agreement between these different versions of
CO2SYS was demonstrated by Orr et al. (2015).  More recently,
CO2SYS-MATLAB has been modified to include uncertainty propagation
(Orr et al., 2018): the main routine CO2SYS.m was altered slightly,
while two new routines were added (errors.m and derivnum.m)

If you discover inconsistencies or have a more general bug report for
CO2SYS.m, please notify S. van Heuven (svheuven at gmail.com), Denis
Pierrot (Denis.Pierrot at noaa.gov), or Alex Kozyr (kozyr at
ornl.gov). For any concerns about the uncertainty propagation routines
(errors.m and derivnum.m), please contact James Orr (james.orr at
lsce.ipsl.fr)

**INSTALLING**

Download the m-files in the src directory (CO2SYS.m, errors.m, and derivnum.m);
you may also wish to download the examples in the examples directory.  Place
these files in a local directory that is in MATLAB's search path, or add the
directory where they are located to MATLAB's search path. The latter can be
done with MATLAB's addpath command, for example

addpath ("my_matlab_directory/my_subdir")

Then run either of the examples in Matlab, or start using the CO2SYS routine
straight away.

**COMPATIBILITY**

Besides their use in MATLAB, the three functions (CO2SYS.m, derivnum.m, and
errors.m) also work well under octave, GNU's MATLAB clone.

**EXAMPLES**

Example MATLAB scripts demonstrating use of CO2SYS can be found in the
examples directory. Using the two new routines is similar, adding only
a few new arguments, e.g., for input uncertainties.  More elaborate
examples are also available in another form in the 'notebooks'
directory. Either click on those files to visualize them (as HTML) or
download those files and use them interactively as jupyter
notebooks. Within MATLAB or octave, you may also use the native help
facility, i.e., by typing "help errors" or "help derivnum" to find out
more.

**REFERENCES**

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2
System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf.  Anal. Cent.,
Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp.,
https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf

Orr, J. C., J.-P. Gattuso, and J.-M. Epitalon (2015) Comparison of ten
packages that compute ocean carbonate chemistry, Biogeosciences, 12,
1483–1510, https://doi.org/10.5194/bg-12-1483-2015 .

Orr, J.C., J.-M. Epitalon, A. G. Dickson, and J.-P. Gattuso (2018) Routine
uncertainty propagation for the marine carbon dioxide system, in prep. for
Mar. Chem., in press, https://doi.org/10.1016/j.marchem.2018.10.006.

van Heuven, S., D. Pierrot, J.W.B. Rae, E. Lewis, and D.W.R. Wallace (2011)
MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b.  Carbon
Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S.
Department of Energy, Oak Ridge, Tennessee. https://doi.org/10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1
=======
# CO2-System-Extd

<a href="https://zenodo.org/badge/latestdoi/198885961"><img src="https://zenodo.org/badge/198885961.svg" alt="DOI"></a> [![View CO2SYSv3 for MATLAB on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/78378-co2sysv3-for-matlab)

## ABOUT

This repository includes software compatible with MATLAB and GNU Octave for calculating marine CO2 system variables (CO2SYS.m), computing partial derivatives of calculated CO2 system variables with respect to inputs (derivnum.m), and propagating uncertainties for CO2 system calculations (errors.m). This software performs similarly to previously released versions of CO2SYS.m (v1: https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/; v2: https://github.com/jamesorr/CO2SYS-MATLAB), and includes the following extended capabilities, additions, and bug fixes (among other minor changes):
 
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

Also included in this repository is a routine to compare CO2SYSv3 to CO2SYSv2 (compare_versions.m), a routine to calculate total concentrations of conservative elements (Na, Mg, Cl, etc.) from CO2SYS output (TOTALS.m), and an example function to run CO2SYSv3 and plot some of the output (example_CO2SYS.m).

## HISTORY

CO2SYS was initially developed by Lewis and Wallace (1998) for MS DOS, later adapted for MS Excel and MATLAB by Pierrot (2006). The code was vectorized, refined, and optimized for computational speed by van Heuven (2011). Options for error propagation were added by Orr et al. (2018). This software builds upon those previous versions.

## INSTALLATION AND USE

Download the files in this repository and place them in a directory that is in the MATLAB search path. Or add the directory where they are located to the search path (https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html).

To perform CO2 system calculations, use CO2SYS.m as directed in the function's help text. All sub-routines that will be called upon are contained within the CO2SYS.m function.

To propagate uncertainties in CO2 system calculations, use errors.m as directed in the function's help text. The errors.m routine will call upon CO2SYS.m and derivnum.m. As such, this version of errors.m is compatible only with CO2SYSv3.

To compute partial derivatives of calculated CO2 system variables with respect to input parameters, use derivnum.m as directed in the function's help text. The derivnum.m routine will call upon CO2SYS.m. As such, this version of derivnum.m is compatible only with CO2SYSv3.

## CITATION

The full citation for CO2SYSv3 (Sharp et al., 2023) is given below. Cite this version if using CO2SYSv3 for CO2 system calculations or propagating errors in CO2 system calculations using the extended errors.m or derivnum.m routines provided here.

If using any CO2SYS program for CO2 system calculations, cite also the original CO2SYS DOS work of Lewis and Wallace (1998).

If using the CO2SYS MATLAB program for CO2 system calculations, cite also the work of van Heuven et al. (2011).

If using the derivnum.m and/or errors.m programs for CO2 system error propagations, cite also the work of Orr et al. (2018).

## REFERENCES

Bach, L. T. (2015). Reconsidering the role of carbonate ion concentration in calcification by marine organisms. Biogeosciences 12(16), 4939–4951.

Dillon, W. D. N., Dillingham, P. W., Currie, K. I., McGraw, C. M., 2020. Inclusion of uncertainty in the calcium-salinity relationship improves estimates of ocean acidification monitoring data quality. Marine Chemistry 226, 103872.

Humphreys, M.P., Lewis, E.R., Sharp, J.D., Pierrot, D. (2022). PyCO2SYS: marine carbonate system calculations in Python. Geoscientific Model Development 15, 15-43..

Lewis, E., Wallace, D. W. R., 1998. Program Developed for CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Munhoven, G., Mathematics of the total alkalinity–pH equation – pathway to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. Geoscientific Model Development 6, 1367–1388, 2013

Orr, J.C., Epitalon, J.-M., Dickson, A. G., Gattuso, J.-P., 2018. Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.

Pierrot, D. E. Lewis,and D. W. R. Wallace. 2006. MS Excel Program Developed for CO2 System Calculations. ORNL/CDIAC-105a. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tennessee. 

Sharp, J.D., Pierrot, D., Humphreys, M.P., Epitalon, J.-M., Orr, J.C., Lewis, E.R., Wallace, D.W.R. (2023, Jan. 19). CO2SYSv3 for MATLAB (Version v3.2.1). Zenodo. http://doi.org/10.5281/zenodo.3950562

Sulpis, O., Lauvset, S. K., and Hagens, M., 2020. Current estimates of K1* and K2* appear inconsistent with measured CO2 system parameters in cold oceanic regions. Ocean Science Discussions, 1-27.

van Heuven, S., Pierrot, D., Rae, J.W.B., Lewis, E., Wallace, D.W.R., 2011. MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.
