HISTORY:
Original version for DOS was written by Lewis and Wallace. This was converted to MATLAB by Denis Pierrot 
at CIMAS, University of Miami, Miami, Florida. Vectorization, internal refinements and speed improvements 
were added by Steven van Heuven, University of Groningen, The Netherlands. Although functionality has 
been added, the output of the function has not changed. All versions of co2sys available at CDIAC 
(DOS, Excel for WINDOWS, Excel for MAC, MATLAB) should yield (near-) identical results when supplied 
with identical input. If you discover that they don't or you have a more general bug report, please 
email me, Denis Pierrot or Alex Kozyr (svheuven@gmail.com, Denis.Pierrot@noaa.gov, kozyra@ornl.gov). 

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

CO2SYS was initially developed by Lewis and Wallace (1998) for MS DOS, later adapted for MS Excel and MATLAB by Pierrot (2006). The code was vectorized, refined, and optimized for computational speed by van Heuven (2011). Options for error propagation were added by Orr et al. (2018). This software builds upon those previous versions.


**CO2SYS-MATLAB versions**

- 1.1   (Sept 2011): van Heuven et al. (2011) 
- 2.0   (20 Dec 2016): includes uncertainty propagation
- 2.0.1 (11 Oct 2017): supports TEOS-10 standards (conservattive temperature, absolute salinity)
- 2.0.2 (17 Oct 2017): Octave enhancements changed to be MATLAB compatible
- 2.0.3 (4 Jun 2018): examples added as Jupyter notebooks
- 2.0.4 (10 Nov 2018): defaults for standard uncertainties in constants (epK vector and eBt) made consistent with Orr et al. (2018), i.e., final published version 
- 2.0.5 (23 Nov 2018): fixed bug in eBt propagation to deriv array (thanks A. Cochon)
- 2.1   (29 Jun 2020): fixed bug in derivnum affecting OUT results (linked to TEMPOUT); masked derivs of input vars in derivnum