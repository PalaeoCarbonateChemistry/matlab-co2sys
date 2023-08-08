[![Run small verification tests](https://github.com/MUADh3i9yL/matlab-co2sys-traditional/actions/workflows/small_tests.yml/badge.svg)](https://github.com/MUADh3i9yL/matlab-co2sys-traditional/actions/workflows/small_tests.yml)
![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fgist.githubusercontent.com%2Frossidae%2Fc68eb447f90f281a543bca7ab1d7a56a%2Fraw%2F9ac172abdacff106379815d86e31720f80dadf57%2Fco2sys-performance-metrics.json&query=%24.iterative&suffix=%20%CE%BCsec%2Fpoint&label=Iterative%20speed&color=%23963023)
![Dynamic JSON Badge](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fgist.githubusercontent.com%2Frossidae%2Fc68eb447f90f281a543bca7ab1d7a56a%2Fraw%2F9ac172abdacff106379815d86e31720f80dadf57%2Fco2sys-performance-metrics.json&query=%24.array&suffix=%20%CE%BCsec%2Fpoint&label=Vectorised%20speed)

# About CO2SYS
This is a Matlab version of CO2SYS (which was originally written for DOS). CO2SYS calculates and returns a the state of the carbonate system of solutions. Give two carbonate system properties and input to calculate apparent equilibrium constants, CO2SYS will return a matrix of results describing subsidiary carbonate system properties.

This is the traditional version of CO2SYS - providing a drag and drop replacement to CO2SYSv3.2, but with improvements to the code quality and speed of execution. CO2SYS accepts arrays of input for vectorised calculation.

# Citation
Whiteford _et al._, 2024, Updating CO2SYS

# Installation
Install from the Matlab file exchange: xxx

Or download a release from GitHub: xxx

Or use Git to install a copy:
```
git init
git remote add origin https://github.com/MUADh3i9yL/matlab-co2sys-traditional
git pull origin main
```

To check installation was successful, run the examples:
```
cd examples
run_examples
```
Ensure that the primary CO2SYS file (in the main folder) is on your path when you run this.

# Examples
Example Matlab scripts demonstrating use of CO2SYS can be found in the examples directory.
Workbooks giving interactive example are found in the Matlab live documents in the same directory.

# Compatibility
Besides their use in Matlab, the three functions (CO2SYS.m, derivnum.m, and errors.m) also work well under octave, GNU's Matlab clone.

# References
Lewis, E., Wallace, D. W. R., 1998. Program Developed for CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Pierrot, D. E. Lewis,and D. W. R. Wallace, 2006, MS Excel Program Developed for CO2 System Calculations. ORNL/CDIAC-105a. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tennessee. 

van Heuven, S., Pierrot, D., Rae, J.W.B., Lewis, E., Wallace, D.W.R., 2011. Matlab Program Developed for CO2 System Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Munhoven, G., 2013, Mathematics of the total alkalinity–pH equation – pathway to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. Geoscientific Model Development 6, 1367–1388

Orr, J. C., J.-P. Gattuso, and J.-M. Epitalon, 2015, Comparison of ten packages that compute ocean carbonate chemistry, Biogeosciences, 12, 1483–1510, https://doi.org/10.5194/bg-12-1483-2015.


Orr, J.C., Epitalon, J.-M., Dickson, A. G., Gattuso, J.-P., 2018, Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.

Humphreys, M.P., Lewis, E.R., Sharp, J.D., Pierrot, D., 2022, PyCO2SYS: marine carbonate system calculations in Python. Geoscientific Model Development 15, 15-43.

Sharp, J.D., Pierrot, D., Humphreys, M.P., Epitalon, J.-M., Orr, J.C., Lewis, E.R., Wallace, D.W.R. (2023, Jan. 19). CO2SYSv3 for Matlab (Version v3.2.1). Zenodo. http://doi.org/10.5281/zenodo.3950562

