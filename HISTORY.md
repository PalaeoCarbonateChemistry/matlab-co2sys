# History
The original version of CO2SYS was written for DOS by Lewis and Wallace. This was converted to MATLAB by Denis Pierrot at CIMAS, University of Miami, Miami, Florida. Vectorization, internal refinements and speed improvements were added by Steven van Heuven, University of Groningen, The Netherlands. Although functionality has been added, the output of the function has not changed.

All versions of co2sys available at CDIAC (DOS, Excel for WINDOWS, Excel for MAC, MATLAB) should yield (near-) identical results when supplied with identical input. Indeed, close agreement between these different versions of CO2SYS was demonstrated by Orr et al. (2015).  More recently, CO2SYS-MATLAB has been modified to include uncertainty propagation (Orr et al., 2018): the main routine CO2SYS.m was altered slightly,while two new routines were added (errors.m and derivnum.m)

# Matlab-CO2SYS versions
- 3.2 (January 2023): included ammonia and sulphide, added new options for K calculation, and fixed various bugs
- 2.1   (29 Jun 2020): fixed bug in derivnum affecting OUT results (linked to TEMPOUT); masked derivs of input vars in derivnum
- 2.0.5 (23 Nov 2018): fixed bug in eBt propagation to deriv array (thanks A. Cochon)
- 2.0.4 (10 Nov 2018): defaults for standard uncertainties in constants (epK vector and eBt) made consistent with Orr et al. (2018), i.e., final published version 
- 2.0.3 (4 Jun 2018): examples added as Jupyter notebooks
- 2.0.2 (17 Oct 2017): Octave enhancements changed to be MATLAB compatible
- 2.0.1 (11 Oct 2017): supports TEOS-10 standards (conservattive temperature, absolute salinity)
- 2.0   (20 Dec 2016): includes uncertainty propagation
- 1.1   (Sept 2011): van Heuven et al. (2011) 





