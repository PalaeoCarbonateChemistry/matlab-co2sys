[![Run small verification tests](https://github.com/MUADh3i9yL/matlab-co2sys-traditional/actions/workflows/small_tests.yml/badge.svg)](https://github.com/MUADh3i9yL/matlab-co2sys-traditional/actions/workflows/small_tests.yml)

![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fgist.githubusercontent.com%2Frossidae%2Fc68eb447f90f281a543bca7ab1d7a56a%2Fraw%2F9ac172abdacff106379815d86e31720f80dadf57%2Fco2sys-performance-metrics.json&query=%24.iterative&suffix=%20%CE%BCsec%2Fpoint&label=Iterative%20speed&color=%23963023)


![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fgist.githubusercontent.com%2Frossidae%2Fc68eb447f90f281a543bca7ab1d7a56a%2Fraw%2F9ac172abdacff106379815d86e31720f80dadf57%2Fco2sys-performance-metrics.json&query=%24.array&suffix=%20%CE%BCsec%2Fpoint&label=Vectorised%20speed)

# About CO2SYS
This is a MATLAB version of CO2SYS (which was originally written for DOS). CO2SYS calculates and returns a the state of the carbonate system of solutions. Give two carbonate system properties and input to calculate apparent equilibrium constants, CO2SYS will return a matrix of results describing subsidiary carbonate system properties.

This is the traditional version of CO2SYS - providing a drag and drop replacement to CO2SYSv3.0, but with improvements to the code quality and speed of execution. CO2SYS accepts arrays of input for vectorised calculation.

# Citation
- If you use any CO2SYS related software, please cite the original work by Lewis and Wallace (1998).
- If you use CO2SYS.m, please cite van Heuven et al (2011).
- If you use errors.m or derivnum.m, please cite Orr et al. (2018).

van Heuven, S., D. Pierrot, J.W.B. Rae, E. Lewis, and D.W.R. Wallace. 2011.
MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b.
Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S.
Department of Energy, Oak Ridge, Tennessee. doi: 10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1

# Installation
Download the m-files in the src directory (CO2SYS.m, errors.m, and derivnum.m);
you may also wish to download the examples in the examples directory.  Place
these files in a local directory that is in MATLAB's search path, or add the
directory where they are located to MATLAB's search path. The latter can be
done with MATLAB's addpath command, for example

addpath ("my_matlab_directory/my_subdir")

Then run either of the examples in Matlab, or start using the CO2SYS routine
straight away.
Download the m-file "CO2SYS.m" and, optionally, the two examples. Place the file(s) in a location that 
Matlab can see (or add the location of the file(s) to Matlab's search path).
Run either of the examples in Matlab, or start using the main routine straight away.

# Examples
Example MATLAB scripts demonstrating use of CO2SYS can be found in the
examples directory. Using the two new routines is similar, adding only
a few new arguments, e.g., for input uncertainties.  More elaborate
examples are also available in another form in the 'notebooks'
directory. Either click on those files to visualize them (as HTML) or
download those files and use them interactively as jupyter
notebooks. Within MATLAB or octave, you may also use the native help
facility, i.e., by typing "help errors" or "help derivnum" to find out
more.

# Compatibility
Besides their use in MATLAB, the three functions (CO2SYS.m, derivnum.m, and
errors.m) also work well under octave, GNU's MATLAB clone.

# References
Bach, L. T. (2015). Reconsidering the role of carbonate ion concentration in calcification by marine organisms. Biogeosciences 12(16), 4939–4951.

Dillon, W. D. N., Dillingham, P. W., Currie, K. I., McGraw, C. M., 2020. Inclusion of uncertainty in the calcium-salinity relationship improves estimates of ocean acidification monitoring data quality. Marine Chemistry 226, 103872.

Humphreys, M.P., Lewis, E.R., Sharp, J.D., Pierrot, D. (2022). PyCO2SYS: marine carbonate system calculations in Python. Geoscientific Model Development 15, 15-43..

Lewis, E., Wallace, D. W. R., 1998. Program Developed for CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Munhoven, G., Mathematics of the total alkalinity–pH equation – pathway to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. Geoscientific Model Development 6, 1367–1388, 2013

Orr, J. C., J.-P. Gattuso, and J.-M. Epitalon (2015) Comparison of ten
packages that compute ocean carbonate chemistry, Biogeosciences, 12,
1483–1510, https://doi.org/10.5194/bg-12-1483-2015 .

Orr, J.C., Epitalon, J.-M., Dickson, A. G., Gattuso, J.-P., 2018. Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.

Pierrot, D. E. Lewis,and D. W. R. Wallace. 2006. MS Excel Program Developed for CO2 System Calculations. ORNL/CDIAC-105a. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tennessee. 

Sharp, J.D., Pierrot, D., Humphreys, M.P., Epitalon, J.-M., Orr, J.C., Lewis, E.R., Wallace, D.W.R. (2023, Jan. 19). CO2SYSv3 for MATLAB (Version v3.2.1). Zenodo. http://doi.org/10.5281/zenodo.3950562

Sulpis, O., Lauvset, S. K., and Hagens, M., 2020. Current estimates of K1* and K2* appear inconsistent with measured CO2 system parameters in cold oceanic regions. Ocean Science Discussions, 1-27.

van Heuven, S., Pierrot, D., Rae, J.W.B., Lewis, E., Wallace, D.W.R., 2011. MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

