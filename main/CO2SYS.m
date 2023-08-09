function [DATA,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,varargin)
%**************************************************************************
%
% Current: CO2SYS.m v3.1.2   (Jan  2023: https://github.com/jonathansharp/CO2-System-Extd)
%          CO2SYS.m v2       (Dec  2016: https://github.com/jamesorr/CO2SYS-MATLAB)
%          CO2SYS.m v1       (Sept 2011: https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/)
%
% CO2SYS is a MATLAB-version of the original CO2SYS for DOS. 
% CO2SYS calculates and returns the state of the carbonate system of 
%    oceanographic water samples, if supplied with enough input.
%
% Please note that, besides certain extended capabilities and minor
%    corrections, this software is intended to be nearly identical to the
%    DOS and Excel versions that have been released previously, meaning
%    that results obtained should be very nearly identical for identical
%    input.
% Additionally, several of the dissociation constants K1 and K2 that have 
%    been published since the original DOS version was written are implemented.
%    For a complete list of changes since version 1.0, see below.
%
% For much more info please have a look at:
%    Lewis, E., and D. W. R. Wallace. 1998. Program Developed for
%    CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information
%    Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy,
%    Oak Ridge, Tennessee. http://cdiac.ornl.gov/oceans/co2rprt.html
%
%**************************************************************************
%
%  **** SYNTAX:
%  [RESULT,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%        ...SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,...
%        ...K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,varargin)
% 
%  **** SYNTAX EXAMPLES:
%  [Result]                     = CO2SYS(2400,2200,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [Result,Headers]             = CO2SYS(2400,   8,1,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [Result,Headers,Niceheaders] = CO2SYS( 500,   8,5,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1,'co2_press',1)
%  [A]                          = CO2SYS(2400,2000:10:2400,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [A]                          = CO2SYS(2400,2200,1,2,0:1:35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [A]                          = CO2SYS(2400,2200,1,2,35,0,25,0:100:4200,0,15,1,0,0,1,4,1,1,1,'co2_press',1)
%  
%  **** APPLICATION EXAMPLE (copy and paste this into command window):
%  tmps=0:40; sals=0:40; [X,Y]=meshgrid(tmps,sals);
%  A = CO2SYS(2300,2100,1,2,Y(:),X(:),nan,0,nan,1,1,0,0,1,9,1,1,1);
%  Z=nan(size(X)); Z(:)=A(:,4); figure; contourf(X,Y,Z,20); caxis([0 1200]); colorbar;
%  ylabel('Salinity [psu]'); xlabel('Temperature [degC]'); title('Dependence of pCO2 [uatm] on T and S')
% 
%**************************************************************************
%
% INPUT:
%
%   PAR1  (some unit) : scalar or vector of size n
%   PAR2  (some unit) : scalar or vector of size n
%   PAR1TYPE       () : scalar or vector of size n (Footnote *1 below)
%   PAR2TYPE       () : scalar or vector of size n (Footnote *1 below)
%   SAL            () : scalar or vector of size n
%   TEMPIN  (degr. C) : scalar or vector of size n 
%   TEMPOUT (degr. C) : scalar or vector of size n 
%   PRESIN     (dbar) : scalar or vector of size n 
%   PRESOUT    (dbar) : scalar or vector of size n
%   SI    (umol/kgSW) : scalar or vector of size n
%   PO4   (umol/kgSW) : scalar or vector of size n
%   NH4   (umol/kgSW) : scalar or vector of size n
%   H2S   (umol/kgSW) : scalar or vector of size n
%   pHSCALEIN         : scalar or vector of size n (Footnote *2 below)
%   K1K2CONSTANTS     : scalar or vector of size n (Footnote *3 below)
%   KSO4CONSTANT      : scalar or vector of size n (Footnote *4 below)
%   KFCONSTANT        : scalar or vector of size n (Footnote *5 below)
%   BORON             : scalar or vector of size n (Footnote *6 below)
%
% OPTIONAL INPUT:
%   'co2_press',X     : string,scalar pair
%               X = 0 : K0 and FugFac are not corrected for in situ pressure
%               X = 1 : K0 and FugFac are corrected for in situ pressure
%
%  (*1) Each element must be an integer, 
%      indicating that PAR1 (or PAR2) is of type: 
%  1 = TA
%  2 = DIC
%  3 = pH
%  4 = pCO2
%  5 = fCO2
%  6 = HCO3
%  7 = CO3
%  8 = CO2
% 
%  (*2) Each element must be an integer, 
%       indicating that the pH-input (PAR1 or PAR2, if any) is at:
%  1 = Total scale
%  2 = Seawater scale
%  3 = Free scale
%  4 = NBS scale
% 
%  (*3) Each element must be an integer, 
%       indicating the K1 and K2 dissociation constants that are to be used:
%   1 = Roy, 1993                                           T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson                                     T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO   T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO   T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO   T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)                   T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng	(i.e., originam Mehrbach but without XXX)   T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)    T:    0-50  S:     0. 
%   9 = Cai and Wang, 1998                                  T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000                                  T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002                     T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002                                 T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006                                 T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero, 2010                                       T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  15 = Waters, Millero, & Woosley, 2014                    T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  16 = Sulpis et al, 2020                                  T: -1.7-32  S: 31-38. Total scale. Field measurements.
%  17 = Schockman & Byrne, 2021                             T:   15-35  S: 19-41. Total scale. Real seawater.
% 
%  (*4) Each element must be an integer that
%       indicates the KSO4 dissociation constant that is to be used:
%  1 = KSO4 of Dickson   (PREFERRED)
%  2 = KSO4 of Khoo
%  3 = KSO4 of Waters and Millero
%
%  (*5) Each element must be an integer that 
%       indicates the KHF dissociation constant that is to be used:
%  1 = KF of Dickson & Riley 1979  
%  2 = KF of Perez & Fraga, 1987  (PREFERRED)
%
%  (*6) Each element must be an integer that 
%       indicates the the formulation of the borate-to-salinity ratio to be used:
%  1 = TB of Uppstrom 1979
%  2 = TB of Lee 2010  (PREFERRED)
%
%**************************************************************************
%
% OUTPUT: *  an array containing the following parameter values (one row per sample):
%         *  a cell-array containing crudely formatted headers
%         *  a cell-array containing nicely formatted headers
%
%   POS   PARAMETER            UNIT
%
%    01 - TAlk                 (umol/kgSW)
%    02 - TCO2                 (umol/kgSW)
%    03 - pHin                 ()
%    04 - pCO2 input           (uatm)
%    05 - fCO2 input           (uatm)
%    06 - HCO3 input           (umol/kgSW)
%    07 - CO3 input            (umol/kgSW)
%    08 - CO2 input            (umol/kgSW)
%    09 - BAlk input           (umol/kgSW)
%    10 - OH input             (umol/kgSW)
%    11 - PAlk input           (umol/kgSW)
%    12 - SiAlk input          (umol/kgSW)
%    13 - AmmAlk input         (umol/kgSW)
%    14 - HSAlk input          (umol/kgSW)
%    15 - Hfree input          (umol/kgSW)
%    16 - RevelleFactor input  ()
%    17 - OmegaCa input        ()
%    18 - OmegaAr input        ()
%    19 - xCO2 input           (ppm)
%    20 - SIR input            ()
%    21 - pH output            ()
%    22 - pCO2 output          (uatm)
%    23 - fCO2 output          (uatm)
%    24 - HCO3 output          (umol/kgSW)
%    25 - CO3 output           (umol/kgSW)
%    26 - CO2 output           (umol/kgSW)
%    27 - BAlk output          (umol/kgSW)
%    28 - OH output            (umol/kgSW)
%    29 - PAlk output          (umol/kgSW)
%    30 - SiAlk output         (umol/kgSW)
%    31 - AmmAlk output        (umol/kgSW)
%    32 - HSAlk output         (umol/kgSW)
%    33 - Hfree output         (umol/kgSW)
%    34 - RevelleFactor output ()
%    35 - OmegaCa output       ()
%    36 - OmegaAr output       ()
%    37 - xCO2 output          (ppm)
%    38 - SIR output            ()
%    39 - pH input (Total)     ()          
%    40 - pH input (SWS)       ()          
%    41 - pH input (Free)      ()          
%    42 - pH input (NBS)       ()          
%    43 - pH output (Total)    ()          
%    44 - pH output (SWS)      ()          
%    45 - pH output (Free)     ()          
%    46 - pH output (NBS)      () 
%    47 - TEMP input           (deg C)     ***    
%    48 - TEMPOUT              (deg C)     ***
%    49 - PRES input           (dbar or m) ***
%    50 - PRESOUT              (dbar or m) ***
%    51 - PAR1TYPE             (integer)   ***
%    52 - PAR2TYPE             (integer)   ***
%    53 - K1K2CONSTANTS        (integer)   ***
%    54 - KSO4CONSTANT         (integer)   *** 
%    55 - KFCONSTANT           (integer)   *** 
%    56 - BORON                (integer)   *** 
%    57 - pHSCALE of input     (integer)   ***
%    58 - SAL                  (psu)       ***
%    59 - PO4                  (umol/kgSW) ***
%    60 - SI                   (umol/kgSW) ***
%    61 - NH4                  (umol/kgSW) ***
%    62 - H2S                  (umol/kgSW) ***
%    63 - K0  input            ()          
%    64 - K1  input            ()          
%    65 - K2  input            ()          
%    66 - pK1 input            ()          
%    67 - pK2 input            ()          
%    68 - KW  input            ()          
%    69 - KB  input            ()          
%    70 - KF  input            ()          
%    71 - KS  input            ()          
%    72 - KP1 input            ()          
%    73 - KP2 input            ()          
%    74 - KP3 input            ()          
%    75 - KSi input            ()              
%    76 - KNH4input            ()    
%    77 - KH2Sinput            ()    
%    78 - K0  output           ()          
%    79 - K1  output           ()          
%    80 - K2  output           ()          
%    81 - pK1 output           ()          
%    82 - pK2 output           ()          
%    83 - KW  output           ()          
%    84 - KB  output           ()          
%    85 - KF  output           ()          
%    86 - KS  output           ()          
%    87 - KP1 output           ()          
%    88 - KP2 output           ()          
%    89 - KP3 output           ()          
%    90 - KSi output           ()              
%    91 - KNH4output           ()    
%    92 - KH2Soutput           () 
%    93 - TB                   (umol/kgSW) 
%    94 - TF                   (umol/kgSW) 
%    95 - TS                   (umol/kgSW) 
%    96 - TP                   (umol/kgSW) 
%    97 - TSi                  (umol/kgSW)
%    98 - TNH4                 (umol/kgSW) 
%    99 - TH2S                 (umol/kgSW)
%
%    *** SIMPLY RESTATES THE INPUT BY USER 
%
% In all the above, the terms "input" and "output" may be understood
%    to refer to the 2 scenarios for which CO2SYS performs calculations, 
%    each defined by its own combination of temperature and pressure.
%    For instance, one may use CO2SYS to calculate, from measured DIC and
%    TAlk, the pH that that sample will have in the lab (e.g., T=25 degC, P=0
%    dbar), and what the in situ pH would have been (e.g., at T=1 degC, P=4500).
%    A = CO2SYS(2400,2200,1,2,35,25,1,0,4200,1,1,0,0,1,4,1,1,1)
%    pH_lab = A(3);  % 7.8429
%    pH_sea = A(20); % 8.0503
% 
%**************************************************************************
%
% **** Changes since 3.1 by JD Sharp.
%   - rigorous validation performed against PyCO2SYS
%     (https://github.com/mvdh7/PyCO2SYS)
%   - initial pH estimates obtained via the approach of Munhoven (2013)
%   - correction to solution for free scale pH within iterative pH solvers
%   - correction to uncertainty calculation for parameters at output conditions
%   - consitency implemented for [CO2(aq)] calculations
%   - substrate-inhibitor ratio (SIR; Bach, 2015) included as an output argument
%   - input uncertainty in [CO2], [HCO3], and [CO3] should now be in mol/kg
%   - option added for pressure corrections to K0 and fugacity factor
%
% **** Changes since 3.0 by JD Sharp.
%   - added KSO4 of Waters and Millero (2013)
%   - added K1 and K2 of Sulpis et al. (2020)
%   - added K1 and K2 of Schockman and Byrne (2021)
%
% **** Changes since 3.0 by JD Sharp based on code from D Pierrot.
%   - changed code to set pH values that don't converge to NaN. All	
%     subsequent calculated values also set to NaN.
%   - modified input function to separate KHSO4 and TB choices
%   - added KHF of Perez & Fraga as choice for HF dissociation constant
%   - modified output to reflect all changes mentioned above
%
% **** Changes since 3.0 by MP Humphreys.
%   - include Peng correction for Icase 16 and 17.
%   - fix Icase typo for CO2-HCO3 input pair.
%   - make corrections to (F) indexing in a few places.
%
% **** Changes since 2.1 (uploaded to GitHub Jul 2019) by JD Sharp
%	- now allows for input of NH4+ and H2S concentrations
%
% **** Additional changes since 2.0 by JD Sharp
%	- now allows for HCO3, CO3, and CO2 as input parameters for calculation and
%     for error propagation
%
% **** Changes since 2.0
%	- slight changes to allow error propagation
%	- new option to choose K1 & K2 from Waters et al. (2014): fixes inconsistencies with Millero (2010) identified by Orr et al. (2015)
%
% **** Changes since 1.01 (uploaded to CDIAC at June 11th, 2009):
% - Function cleans up its global variables when done (if you lose variables, this may be the cause -- see around line 570)
% - Added the outputting of K values
% - Implementation of constants of Cai and Wang, 1998
% - Implementation of constants of Lueker et al., 2000
% - Implementation of constants of Mojica-Prieto and Millero, 2002
% - Implementation of constants of Millero et al., 2002 (only their eqs. 19, 20, no TCO2 dependency)
% - Implementation of constants of Millero et al., 2006
% - Implementation of constants of Millero et al., 2010
% - Properly listed Sal and Temp limits for the available constants
% - added switch for using the new Lee et al., (2010) formulation of Total Borate
% - Minor corrections to the GEOSECS constants (gave NaN for some output in earlier version)
% - Fixed decimal point error on [H+] (did not get converted to umol/kgSW from mol/kgSW).
% - Changed 'Hfreein' to 'Hfreeout' in the 'NICEHEADERS'-output (typo)
%
% **** Changes since 1.00 (uploaded to CDIAC at May 29th, 2009):
% - added a note explaining that all known bugs were removed before release of 1.00
%
%**************************************************************************
%
% CO2SYS originally by Lewis and Wallace 1998
%
% Converted to MATLAB by Denis Pierrot at
% CIMAS, University of Miami, Miami, Florida
%
% Vectorization, internal refinements and speed improvements by
% Steven van Heuven, University of Groningen, The Netherlands.
% Questions, bug reports et cetera: svheuven@gmail.com
%
% Modifications for error propagation by JM Epitalon
%
% Extension to include input of CO2, HCO3, CO3, NH4, and H2S by
% Jonathan Sharp, University of South Florida.
%
% Modification to set pH values that do not converge to NaN, separate
% KHSO4 and TB, and to add the KHF of Perez & Fraga by Denis Pierrot,
% implemented in this version by Jonathan Sharp, University of Washington
%
% Bug fixes by Matthew Humphreys, NIOZ Texel, the Netherlands.
%
% Additional modifications for consistency with PyCO2SYS and other added
% options and input/output arguments by Jonathan Sharp, University of
% Washington
%
%**************************************************************************


%**************************************************************************
% NOTHING BELOW THIS SHOULD REQUIRE EDITING BY USER!
%**************************************************************************


% Declare global variables
global pHScale WhichKs WhoseKSO4 WhoseKF WhoseTB Pbar
global Sal sqrSal TempK logTempK TempCi TempCo Pdbari Pdbaro;
global FugFac VPFac PengCorrection ntps RGasConstant;
global fH RT;
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S CAL F;

% Added by JM Epitalon
% For computing derivative with respect to Ks, one has to call CO2sys with a perturbed K
% Requested perturbation is passed through the following global variables
global PertK    % Id of perturbed K
global Perturb  % perturbation

% Input conditioning

% set default for optional input argument
global p_opt
p_opt = 0;
% parse optional input argument
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'co2_press')
        p_opt = varargin{i+1};
    end
end

% Determine lengths of input vectors
veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
            length(PAR2TYPE) length(SAL) length(TEMPIN)...
            length(TEMPOUT) length(PRESIN) length(PRESOUT)...
            length(SI) length(PO4) length(NH4) length(H2S)...
            length(pHSCALEIN) length(K1K2CONSTANTS) length(KSO4CONSTANT)...
	        length(KFCONSTANT) length(BORON)];

if length(unique(veclengths))>2
	disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
end

% Make column vectors of all input vectors
PAR1         =PAR1         (:);
PAR2         =PAR2         (:);
PAR1TYPE     =PAR1TYPE     (:);
PAR2TYPE     =PAR2TYPE     (:);
SAL          =SAL          (:);
TEMPIN       =TEMPIN       (:);
TEMPOUT      =TEMPOUT      (:);
PRESIN       =PRESIN       (:);
PRESOUT      =PRESOUT      (:);
SI           =SI           (:);
PO4          =PO4          (:);
NH4          =NH4          (:);
H2S          =H2S          (:);
pHSCALEIN    =pHSCALEIN    (:);
K1K2CONSTANTS=K1K2CONSTANTS(:);
KSO4CONSTANT =KSO4CONSTANT (:);
KFCONSTANT   =KFCONSTANT   (:);
BORON        =BORON        (:);

% Find the longest column vector:
ntps = max(veclengths);

% Populate column vectors
PAR1(1:ntps,1)          = PAR1(:)          ;
PAR2(1:ntps,1)          = PAR2(:)          ;
PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
SAL(1:ntps,1)           = SAL(:)           ;
TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
PRESIN(1:ntps,1)        = PRESIN(:)        ;
PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
SI(1:ntps,1)            = SI(:)            ;
PO4(1:ntps,1)           = PO4(:)           ;
NH4(1:ntps,1)           = NH4(:)           ;
H2S(1:ntps,1)           = H2S(:)           ;
pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
KSO4CONSTANT(1:ntps,1)  = KSO4CONSTANT(:)  ;
KFCONSTANT(1:ntps,1)    = KFCONSTANT(:)    ;
BORON(1:ntps,1)         = BORON(:)         ;

% Assign input to the 'historical' variable names.
pHScale      = pHSCALEIN;
WhichKs      = K1K2CONSTANTS;
WhoseKSO4    = KSO4CONSTANT;
WhoseKF      = KFCONSTANT;
WhoseTB      = BORON;
p1           = PAR1TYPE;
p2           = PAR2TYPE;
TempCi       = TEMPIN;
TempCo       = TEMPOUT;
Pdbari       = PRESIN;
Pdbaro       = PRESOUT;
Sal          = SAL;
sqrSal       = sqrt(SAL);
TP           = PO4;
TSi          = SI;
TNH4         = NH4;
TH2S         = H2S;

RGasConstant = 83.14462618; % ml bar-1 K-1 mol-1,
%                             recommended by NIST
%                             https://physics.nist.gov/cgi-bin/cuu/Value?r
%RGasConstant = 83.1451;  % ml bar-1 K-1 mol-1, DOEv2
%                           Compatible w/ CO2SYSv2.0.5
%RGasConstant = 83.14472; % ml bar-1 K-1 mol-1, DOEv3

% Generate empty vectors for...
TA   = nan(ntps,1); % Talk
TC   = nan(ntps,1); % DIC
PH   = nan(ntps,1); % pH
PC   = nan(ntps,1); % pCO2
FC   = nan(ntps,1); % fCO2
HCO3 = nan(ntps,1); % [HCO3]
CO3  = nan(ntps,1); % [CO3]
CO2  = nan(ntps,1); % [CO2*]

% Assign values to empty vectors.
F=(p1==1 & PAR1~=-999);   TA(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==2 & PAR1~=-999);   TC(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==3 & PAR1~=-999);   PH(F)=PAR1(F);
F=(p1==4 & PAR1~=-999);   PC(F)=PAR1(F)/1e6; % Convert from microatm. to atm.
F=(p1==5 & PAR1~=-999);   FC(F)=PAR1(F)/1e6; % Convert from microatm. to atm.
F=(p1==6 & PAR1~=-999); HCO3(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==7 & PAR1~=-999);  CO3(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==8 & PAR1~=-999);  CO2(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==1 & PAR2~=-999);   TA(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==2 & PAR2~=-999);   TC(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==3 & PAR2~=-999);   PH(F)=PAR2(F);
F=(p2==4 & PAR2~=-999);   PC(F)=PAR2(F)/1e6; % Convert from microatm. to atm.
F=(p2==5 & PAR2~=-999);   FC(F)=PAR2(F)/1e6; % Convert from microatm. to atm.
F=(p2==6 & PAR2~=-999); HCO3(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==7 & PAR2~=-999);  CO3(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==8 & PAR2~=-999);  CO2(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg

% Generate the columns holding Si, Phos, Amm, H2S and Sal.
% Pure Water case:
F=(WhichKs==8);
Sal(F) = 0;
% GEOSECS and Pure Water:
F=(WhichKs==8 | WhichKs==6);  
TP(F)  = 0;
TSi(F) = 0;
TNH4(F)  = 0;
TH2S(F)  = 0;
% All other cases
F=~F;                         
TP(F)   = TP(F)./1e6;
TSi(F)  = TSi(F)./1e6;
TNH4(F) = TNH4(F)./1e6;
TH2S(F) = TH2S(F)./1e6;

% The vector 'PengCorrection' is used to modify the value of TA, for those
% cases where WhichKs==7, since PAlk(Peng) = PAlk(Dickson) + TP.
% Thus, PengCorrection is 0 for all cases where WhichKs is not 7
PengCorrection=zeros(ntps,1); F=WhichKs==7; PengCorrection(F)=TP(F);

% Calculate the constants for all samples at input conditions
% The constants calculated for each sample will be on the appropriate pH scale!
calculate_equilibrium_constants(TempCi,Pdbari);

% Added by JM Epitalon
% For computing derivative with respect to Ks, one has to perturb the value of one K
% Requested perturbation is passed through global variables PertK and Perturb
if (~ isempty(PertK))
    switch PertK
        case {'K0'}
            K0 = K0 + Perturb;
        case {'K1'}
            K1 = K1 + Perturb;
        case {'K2'}
            K2 = K2 + Perturb;
        case {'KB'}
            KB = KB + Perturb;
        case {'KW'}
            KW = KW + Perturb;
        case {'BOR'}
            TB = TB + Perturb;
    end
end


% Make sure fCO2 is available for each sample that has pCO2 or CO2.
F = (~isnan(PC) & (p1==4 | p2==4));  FC(F) = PC(F).*FugFac(F);
F = (~isnan(CO2) & (p1==8 | p2==8)); FC(F) = CO2(F)./K0(F);

% Generate vectors for results, and copy the raw input values into them
TAc    = TA;
TCc    = TC;
PHic   = PH;
PCic   = PC;
FCic   = FC;
HCO3ic = HCO3;
CO3ic  = CO3;
CO2ic  = CO2;

% Generate vector describing the combination of input parameters
% So, the valid ones are:
% 12,13,15,16,17,18,23,25,26,27,28,35,36,37,38,56,57,67,68,78
Icase = 10*min(p1,p2) + max(p1,p2);

% Calculate missing values for AT,CT,PH,FC,HCO3,CO3,CO2:
% pCO2 will be calculated later on, routines work with fCO2.
F=Icase==12; % input TA, TC
if any(F)
F=(~isnan(TAc) & ~isnan(TCc) & F);
    PHic(F)                = CalculatepHfromTATC(TAc(F)-PengCorrection(F),TCc(F));
    F=(~isnan(PHic) & F);
    if any(F)
       FCic(F)              = CalculatefCO2fromTCpH(TCc(F), PHic(F));
       [CO3ic(F),HCO3ic(F)] = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==13; % input TA, pH
if any(F)
F=(~isnan(TAc) & ~isnan(PHic) & F);
    TCc(F)                 = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    [CO3ic(F),HCO3ic(F)]   = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==14 | Icase==15 | Icase==18; % input TA, (pCO2 or fCO2 or CO2)
if any(F)
F=(~isnan(TAc) & ~isnan(FCic) & F);
    PHic(F)                = CalculatepHfromTAfCO2(TAc(F)-PengCorrection(F),FCic(F));
    F=(~isnan(PHic) & F);
    if any(F)
       TCc(F)              = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
       [CO3ic(F),HCO3ic(F)]= CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==16; % input TA, HCO3
if any(F)
F=(~isnan(TAc) & ~isnan(HCO3ic) & F);
    PHic(F)                = CalculatepHfromTAHCO3(TAc(F)-PengCorrection(F),HCO3ic(F));  % added Peng correction // MPH
    F=(~isnan(PHic) & F);
    if any(F)
       TCc(F)              = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
       FCic(F)             = CalculatefCO2fromTCpH(TCc(F),PHic(F)); 
       CO3ic(F)            = CalculateCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==17; % input TA, CO3
if any(F)
F=(~isnan(TAc) & ~isnan(CO3ic) & F);
    PHic(F)                 = CalculatepHfromTACO3(TAc(F)-PengCorrection(F),CO3ic(F));  % added Peng correction // MPH
    F=(~isnan(PHic) & F);
    if any(F)
       TCc(F)               = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
       FCic(F)              = CalculatefCO2fromTCpH(TCc(F),PHic(F)); 
       HCO3ic(F)            = CalculateHCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==23; % input TC, pH
if any(F)
F=(~isnan(TCc) & ~isnan(PHic) & F);
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    [CO3ic(F),HCO3ic(F)]    = CalculateCO3HCO3fromTCpH(TCc(F), PHic(F));
end
F=Icase==24 | Icase==25 | Icase==28;  % input TC, (pCO2 or fCO2 or CO2)
if any(F)
F=(~isnan(TCc) & ~isnan(FCic) & F);
    PHic(F)                 = CalculatepHfromTCfCO2(TCc(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    [CO3ic(F),HCO3ic(F)]    = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==26; % input TC, HCO3
if any(F)
F=(~isnan(TCc) & ~isnan(HCO3ic) & F);
    [PHic(F),FCic(F)]       = CalculatepHfCO2fromTCHCO3(TCc(F),HCO3ic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    CO3ic(F)                = CalculateCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==27; % input TC, CO3
if any(F)
F=(~isnan(TCc) & ~isnan(CO3ic) & F);
    [PHic(F),FCic(F)]       = CalculatepHfCO2fromTCCO3(TCc(F),CO3ic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    HCO3ic(F)               = CalculateHCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==34 | Icase==35 | Icase==38; % input pH, (pCO2 or fCO2 or CO2)
if any(F)
F=(~isnan(PHic) & ~isnan(FCic) & F);
    TCc(F)                  = CalculateTCfrompHfCO2(PHic(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    [CO3ic(F),HCO3ic(F)]    = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==36; % input pH, HCO3
if any(F)
F=(~isnan(PHic) & ~isnan(HCO3ic) & F);
    TAc(F)                  = CalculateTAfrompHHCO3(PHic(F),HCO3ic(F)) + PengCorrection(F);
    TCc(F)                  = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    CO3ic(F)                = CalculateCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==37; % input pH, CO3
if any(F)
F=(~isnan(PHic) & ~isnan(CO3ic) & F);
    TAc(F)                  = CalculateTAfrompHCO3(PHic(F),CO3ic(F)) + PengCorrection(F);
    TCc(F)                  = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    HCO3ic(F)               = CalculateHCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==46 | Icase==56 | Icase==68; % input (pCO2 or fCO2 or CO2), HCO3
if any(F)
F=(~isnan(FCic) & ~isnan(HCO3ic) & F);
    PHic(F)                 = CalculatepHfromfCO2HCO3(FCic(F),HCO3ic(F));
    TCc(F)                  = CalculateTCfrompHfCO2(PHic(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    CO3ic(F)                = CalculateCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==47 | Icase==57 | Icase==78; % input (pCO2 or fCO2 or CO2), CO3
if any(F)
F=(~isnan(FCic) & ~isnan(CO3ic) & F);
    PHic(F)                 = CalculatepHfromfCO2CO3(FCic(F),CO3ic(F));
    TCc(F)                  = CalculateTCfrompHfCO2 (PHic(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    HCO3ic(F)               = CalculateHCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==67; % input HCO3, CO3
if any(F)
F=(~isnan(HCO3ic) & ~isnan(CO3ic) & F);
    PHic(F)                 = CalculatepHfromCO3HCO3(CO3ic(F),HCO3ic(F));
    TAc(F)                  = CalculateTAfrompHCO3(PHic(F),CO3ic(F)) + PengCorrection(F);
    TCc(F)                  = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    %CO2ic(F)                = CalculateCO2fromTCpH(TCc(F),PHic(F));
end

% By now, an fCO2 value is available for each sample.
% Generate the associated pCO2 value:
F = (isnan(PCic) & (p1~=4 | p2~=4)); PCic(F)  = FCic(F)./FugFac(F);
% Generate the associated CO2 value:
F = (isnan(CO2ic) & (p1~=8 | p2~=8)); CO2ic(F) = FCic(F).*K0(F);

% Calculate Other Params At Input Conditions:
BAlkinp    = nan(ntps,1); % Generate empty vectors
[OHinp,PAlkinp,SiAlkinp,AmmAlkinp,HSAlkinp,Hfreeinp,HSO4inp,HFinp,...
    Revelleinp,OmegaCainp,OmegaArinp,xCO2dryinp] = deal(BAlkinp);
F=(~isnan(PHic)); % if PHic = NaN, pH calculation was not performed or did not converge
[BAlkinp(F),OHinp(F), PAlkinp(F),SiAlkinp(F),AmmAlkinp(F),...
    HSAlkinp(F), Hfreeinp(F),HSO4inp(F),HFinp(F)] = CalculateAlkParts(PHic(F));
PAlkinp(F)                = PAlkinp(F)+PengCorrection(F);
Revelleinp(F)             = RevelleFactor(TAc(F)-PengCorrection(F), TCc(F));
[OmegaCainp(F),OmegaArinp(F)] = CaSolubility(Sal(F), TempCi(F), Pdbari(F), TCc(F), PHic(F));
xCO2dryinp(~isnan(PCic),1) = PCic(~isnan(PCic),1)./VPFac(~isnan(PCic),1); % ' this assumes pTot = 1 atm
SIRinp = HCO3ic./(Hfreeinp.*1e6);

% % Just for reference, convert pH at input conditions to the other scales
pHicT = nan(ntps,1);
pHicS = nan(ntps,1);
pHicF = nan(ntps,1);
pHicN = nan(ntps,1);
[pHicT(F),pHicS(F),pHicF(F),pHicN(F)]=FindpHOnAllScales(PHic(F));

% Merge the Ks at input into an array. Ks at output will be glued to this later.
KIVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S];

% Calculate the constants for all samples at output conditions
calculate_equilibrium_constants(TempCo,Pdbaro);

% Added by JM Epitalon
% For computing derivative with respect to Ks, one has to perturb the value of one K
% Requested perturbation is passed through global variables PertK and Perturb
if (~ isempty(PertK))
    switch PertK
        case {'K0'}
            K0 = K0 + Perturb;
        case {'K1'}
            K1 = K1 + Perturb;
        case {'K2'}
            K2 = K2 + Perturb;
        case {'KB'}
            KB = KB + Perturb;
        case {'KW'}
            KW = KW + Perturb;
        case {'BOR'}
            TB = TB + Perturb;
    end
end                  

% For output conditions, using conservative TA and TC, calculate pH, fCO2
% and pCO2, HCO3, CO3, and CO2
F=(~isnan(TAc) & ~isnan(TCc)); % i.e., do for all samples that have TA and TC values
PHoc=nan(ntps,1);
[CO3oc,HCO3oc,FCoc] = deal(PHoc);
PHoc(F) = CalculatepHfromTATC(TAc(F)-PengCorrection(F), TCc(F)); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
    FCoc(F) = CalculatefCO2fromTCpH(TCc(F), PHoc(F));
    [CO3oc(F),HCO3oc(F)] = CalculateCO3HCO3fromTCpH(TCc(F),PHoc(F));

% Generate the associated pCO2 value:
PCoc  = FCoc./FugFac;
% Generate the associated CO2 value:
CO2oc = FCoc.*K0;

% Calculate Other Params At Output Conditions:
BAlkout    = nan(ntps,1); % Generate empty vectors
[OHout,PAlkout,SiAlkout,AmmAlkout,HSAlkout,Hfreeout,HSO4out,HFout,...
    Revelleout,OmegaCaout,OmegaArout,xCO2dryout] = deal(BAlkout);
F=(~isnan(PHoc)); % if PHoc = NaN, pH calculation was not performed or did not converge
[BAlkout(F),OHout(F),PAlkout(F),SiAlkout(F),AmmAlkout(F),...
    HSAlkout(F), Hfreeout(F),HSO4out(F),HFout(F)] = CalculateAlkParts(PHoc(F));
PAlkout(F)                 = PAlkout(F)+PengCorrection(F);
Revelleout(F)              = RevelleFactor(TAc(F)-PengCorrection(F), TCc(F));
[OmegaCaout(F),OmegaArout(F)] = CaSolubility(Sal(F), TempCo(F), Pdbaro(F), TCc(F), PHoc(F));
xCO2dryout(~isnan(PCoc),1)    = PCoc(~isnan(PCoc))./VPFac(~isnan(PCoc)); % ' this assumes pTot = 1 atm
SIRout = HCO3oc./(Hfreeout.*1e6);

% Just for reference, convert pH at output conditions to the other scales
pHocT = nan(ntps,1);
pHocS = nan(ntps,1);
pHocF = nan(ntps,1);
pHocN = nan(ntps,1);
[pHocT(F),pHocS(F),pHocF(F),pHocN(F)]=FindpHOnAllScales(PHoc(F));

KOVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S];
TVEC =[TB TF TS TP TSi TNH4 TH2S];

% Saving data in array, 99 columns, as many rows as samples input
DATA=[TAc*1e6         TCc*1e6        PHic           PCic*1e6        FCic*1e6...
      HCO3ic*1e6      CO3ic*1e6      CO2ic*1e6      BAlkinp*1e6     OHinp*1e6...
      PAlkinp*1e6     SiAlkinp*1e6   AmmAlkinp*1e6  HSAlkinp*1e6    Hfreeinp*1e6... %%% Multiplied Hfreeinp *1e6, svh20100827
      Revelleinp      OmegaCainp     OmegaArinp     xCO2dryinp*1e6  SIRinp...
      PHoc            PCoc*1e6       FCoc*1e6       HCO3oc*1e6      CO3oc*1e6...
      CO2oc*1e6       BAlkout*1e6    OHout*1e6      PAlkout*1e6     SiAlkout*1e6...
      AmmAlkout*1e6   HSAlkout*1e6   Hfreeout*1e6   Revelleout      OmegaCaout... %%% Multiplied Hfreeout *1e6, svh20100827
      OmegaArout      xCO2dryout*1e6 SIRout         pHicT           pHicS...
      pHicF           pHicN          pHocT          pHocS           pHocF...
      pHocN           TEMPIN         TEMPOUT        PRESIN          PRESOUT...
      PAR1TYPE        PAR2TYPE       K1K2CONSTANTS  KSO4CONSTANT    KFCONSTANT...
      BORON           pHSCALEIN      SAL            PO4             SI...
      NH4             H2S            KIVEC          KOVEC           TVEC*1e6];
DATA(isnan(DATA))=-999;

HEADERS={'TAlk';'TCO2';'pHin';'pCO2in';'fCO2in';'HCO3in';'CO3in';...
    'CO2in';'BAlkin';'OHin';'PAlkin';'SiAlkin';'AmmAlkin';'HSAlkin';...
    'Hfreein';'RFin';'OmegaCAin';'OmegaARin';'xCO2in';'SIRin';'pHout';...
    'pCO2out';'fCO2out';'HCO3out';'CO3out';'CO2out';'BAlkout';'OHout';...
    'PAlkout';'SiAlkout';'AmmAlkout';'HSAlkout';'Hfreeout';'RFout';'OmegaCAout';...
    'OmegaARout';'xCO2out';'SIRout';'pHinTOTAL';'pHinSWS';'pHinFREE';'pHinNBS';...
    'pHoutTOTAL';'pHoutSWS';'pHoutFREE';'pHoutNBS';'TEMPIN';'TEMPOUT';...
    'PRESIN';'PRESOUT';'PAR1TYPE';'PAR2TYPE';'K1K2CONSTANTS';'KSO4CONSTANT';... KSO4CONSTANTS => KSO4CONSTANT // MPH
    'KFCONSTANT';'BORON';'pHSCALEIN';'SAL';'PO4';'SI';'NH4';'H2S';'K0input';...
    'K1input';'K2input';'pK1input';'pK2input';'KWinput';'KBinput';'KFinput';...
    'KSinput';'KP1input';'KP2input';'KP3input';'KSiinput';'KNH4input';...
    'KH2Sinput';'K0output';'K1output';'K2output';'pK1output';'pK2output';...
    'KWoutput';'KBoutput';'KFoutput';'KSoutput';'KP1output';'KP2output';...
    'KP3output';'KSioutput';'KNH4output';'KH2Soutput';'TB';'TF';'TS';...
    'TP';'TSi';'TNH4';'TH2S'};

NICEHEADERS={...
    '01 - TAlk             (umol/kgSW) ';
    '02 - TCO2             (umol/kgSW) ';
    '03 - pHin             ()          ';
    '04 - pCO2in           (uatm)      ';
    '05 - fCO2in           (uatm)      ';
    '06 - HCO3in           (umol/kgSW) ';
    '07 - CO3in            (umol/kgSW) ';
    '08 - CO2in            (umol/kgSW) ';
    '09 - BAlkin           (umol/kgSW) ';
    '10 - OHin             (umol/kgSW) ';
    '11 - PAlkin           (umol/kgSW) ';
    '12 - SiAlkin          (umol/kgSW) ';
    '13 - AmmAlkin         (umol/kgSW) ';
    '14 - HSAlkin          (umol/kgSW) ';
    '15 - Hfreein          (umol/kgSW) ';
    '16 - RevelleFactorin  ()          ';
    '17 - OmegaCain        ()          ';
    '18 - OmegaArin        ()          ';
    '19 - xCO2in           (ppm)       ';
    '20 - SIRin            ()          ';
    '21 - pHout            ()          ';
    '22 - pCO2out          (uatm)      ';
    '23 - fCO2out          (uatm)      ';
    '24 - HCO3out          (umol/kgSW) ';
    '25 - CO3out           (umol/kgSW) ';
    '26 - CO2out           (umol/kgSW) ';
    '27 - BAlkout          (umol/kgSW) ';
    '28 - OHout            (umol/kgSW) ';
    '29 - PAlkout          (umol/kgSW) ';
    '30 - SiAlkout         (umol/kgSW) ';
    '31 - AmmAlkout        (umol/kgSW) ';
    '32 - HSAlkout         (umol/kgSW) ';
    '33 - Hfreeout         (umol/kgSW) ';
    '34 - RevelleFactorout ()          ';
    '35 - OmegaCaout       ()          ';
    '36 - OmegaArout       ()          ';
    '37 - xCO2out          (ppm)       ';
    '38 - SIRout           ()          ';
    '39 - pHin (Total)     ()          ';
    '40 - pHin (SWS)       ()          ';
    '41 - pHin (Free)      ()          ';
    '42 - pHin (NBS )      ()          ';
    '43 - pHout(Total)     ()          ';
    '44 - pHout(SWS)       ()          ';
    '45 - pHout(Free)      ()          ';
    '46 - pHout(NBS )      ()          ';
    '47 - TEMPIN           (Deg C)     ';    
    '48 - TEMPOUT          (Deg C)     ';
    '49 - PRESIN           (dbar)      ';
    '50 - PRESOUT          (dbar)      ';
    '51 - PAR1TYPE         ()          ';
    '52 - PAR2TYPE         ()          ';
    '53 - K1K2CONSTANTS    ()          ';
    '54 - KSO4CONSTANT     ()          ';
    '55 - KFCONSTANT       ()          ';
    '56 - BORON            ()          ';
    '57 - pHSCALEIN        ()          ';
    '58 - SAL              (umol/kgSW) ';
    '59 - PO4              (umol/kgSW) ';
    '60 - SI               (umol/kgSW) ';
    '61	- NH4	           (umol/kgSW) ';
    '62	- H2S	           (umol/kgSW) ';
    '63 - K0input          ()          ';
    '64 - K1input          ()          ';
    '65 - K2input          ()          ';
    '66 - pK1input         ()          ';
    '67 - pK2input         ()          ';
    '68 - KWinput          ()          ';
    '69 - KBinput          ()          ';
    '70 - KFinput          ()          ';
    '71 - KSinput          ()          ';
    '72 - KP1input         ()          ';
    '73 - KP2input         ()          ';
    '74 - KP3input         ()          ';
    '75 - KSiinput         ()          ';
    '76 - KNH4input        ()          ';
    '77 - KH2Sinput        ()          ';  
    '78 - K0output         ()          ';
    '79 - K1output         ()          ';
    '80 - K2output         ()          ';
    '81 - pK1output        ()          ';
    '82 - pK2output        ()          ';
    '83 - KWoutput         ()          ';
    '84 - KBoutput         ()          ';
    '85 - KFoutput         ()          ';
    '86 - KSoutput         ()          ';
    '87 - KP1output        ()          ';
    '88 - KP2output        ()          ';
    '89 - KP3output        ()          ';
    '90 - KSioutput        ()          ';
    '91 - KNH4output       ()          ';
    '92 - KH2Soutput       ()          ';
    '93 - TB               (umol/kgSW) ';
    '94 - TF               (umol/kgSW) ';
    '95 - TS               (umol/kgSW) ';
    '96 - TP               (umol/kgSW) ';
    '97 - TSi              (umol/kgSW) ';
    '98 - TNH4             (umol/kgSW) ';
    '99 - TH2S             (umol/kgSW) '};

clear global F K2 KP3 Pdbari Sal TS VPFac ntps 
clear global FugFac KB KS Pdbaro T TSi BORON WhichKs pHScale 
clear global K KF KSi KNH4 KH2S PengCorrection TB TempCi WhoseKSO4 WhoseKF WhoseTB sqrSal 
clear global K0 KP1 KW RGasConstant TF TempCo fH 
clear global K1 KP2 Pbar RT TP TempK logTempK
	
end % end main function


%**************************************************************************
% Subroutines:
%**************************************************************************



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core input functions, based on functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculatepHfromTATC(TAi, TCi)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F;
% ' SUB CalculatepHfromTATC, version 04.01, 10-13-96, written by Ernie Lewis
% ' with modifications from Denis Pierrot.
% ' Inputs: TA, TC, K(), T()
% ' Output: pH
%
% ' This calculates pH from TA and TC using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013).
%
% ' Made this to accept vectors. It will continue iterating until all
% ' values in the vector are "abs(deltapH) < pHTol". SVH2007
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input. JDS2020
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl          = sum(F);  % VectorLength
% Find initital pH guess using method of Munhoven (2013)
pHGuess         = CalculatepHfromTATCMunhoven(TAi, TCi);
ln10            = log(10);
pH              = pHGuess;
pHTol           = 0.0001;  % tolerance for iterations end
deltapH(1:vl,1) = pHTol+1;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    Denom     = (H.*H + K1F.*H + K1F.*K2F);
    CAlk      = TCi.*K1F.*(H + 2.*K2F)./Denom;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); % since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree); % since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk  - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % find Slope dTA/dpH;
    % (this is not exact, but keeps all important terms);
    Slope     = ln10.*(TCi.*K1F.*H.*(H.*H + K1F.*K2F + 4.*H.*K2F)./Denom./Denom + BAlk.*H./(KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculatefCO2fromTCpH(TCx, pHx)
global K0 K1 K2 F
% ' SUB CalculatefCO2fromTCpH, version 02.02, 12-13-96, written by Ernie Lewis.
% ' Inputs: TC, pH, K0, K1, K2
% ' Output: fCO2
% ' This calculates fCO2 from TC and pH, using K0, K1, and K2.
H            = 10.^(-pHx);
fCO2x        = TCx.*H.*H./(H.*H + K1(F).*H + K1(F).*K2(F))./K0(F);
varargout{1} = fCO2x;
end % end nested function

function varargout=CalculateTCfromTApH(TAx, pHx)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
% ' SUB CalculateTCfromTApH, version 02.03, 10-10-97, written by Ernie Lewis.
% ' Inputs: TA, pH, K(), T()
% ' Output: TC
% ' This calculates TC from TA and pH.
H         = 10.^(-pHx);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHx); % this converts pH to pHfree no matter the scale
Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree); %' since KF is on the free scale
CAlk      = TAx - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
TCxtemp   = CAlk.*(H.*H + K1F.*H + K1F.*K2F)./(K1F.*(H + 2.*K2F));
varargout{1} = TCxtemp;
end % end nested function

function varargout=CalculatepHfromTAfCO2(TAi, fCO2i)
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculatepHfromTAfCO2, version 04.01, 10-13-97, written by Ernie
% ' Lewis with modifications from Denis Pierrot.
% ' Inputs: TA, fCO2, K0, K(), T()
% ' Output: pH
% ' This calculates pH from TA and fCO2 using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013) and extended by Humphreys et al. (2021).
%
% ' This will continue iterating until all values in the vector are
% ' "abs(deltapH) < pHTol"
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input.
K0F=K0(F);     K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl         = sum(F); % vectorlength
% Find initital pH guess using method of Munhoven (2013)
CO2i       = fCO2i.*K0F; % Convert fCO2 to CO2
pHGuess    = CalculatepHfromTACO2Munhoven(TAi, CO2i);
ln10       = log(10);
pH         = pHGuess;
pHTol      = 0.0001; % tolerance
deltapH = pHTol+pH;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    HCO3      = K0F.*K1F.*fCO2i./H;
    CO3       = K0F.*K1F.*K2F.*fCO2i./(H.*H);
    CAlk      = HCO3 + 2.*CO3;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % '               find Slope dTA/dpH
    % '               (this is not exact, but keeps all important terms):
    Slope     = ln10.*(HCO3 + 4.*CO3 + BAlk.*H./(KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculateTAfromTCpH(TCi, pHi)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculateTAfromTCpH, version 02.02, 10-10-97, written by Ernie Lewis.
% ' Inputs: TC, pH, K(), T()
% ' Output: TA
% ' This calculates TA from TC and pH.
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
H         = 10.^(-pHi);
CAlk      = TCi.*K1F.*(H + 2.*K2F)./(H.*H + K1F.*H + K1F.*K2F);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHi); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
TActemp    = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
varargout{1}=TActemp;
end % end nested function

function varargout=CalculatepHfromTCfCO2(TCi, fCO2i)
global K0 K1 K2 F;
% ' SUB CalculatepHfromTCfCO2, version 02.02, 11-12-96, written by Ernie Lewis.
% ' Inputs: TC, fCO2, K0, K1, K2
% ' Output: pH
% ' This calculates pH from TC and fCO2 using K0, K1, and K2 by solving the
% '       quadratic in H: fCO2.*K0 = TC.*H.*H./(K1.*H + H.*H + K1.*K2).
% ' if there is not a real root, then pH is returned as missingn.
RR = K0(F).*fCO2i./TCi;
%       if RR >= 1
%          varargout{1}= missingn;
%          disp('nein!');return;
%       end
% check after sub to see if pH = missingn.
Discr = (K1(F).*RR).*(K1(F).*RR) + 4.*(1 - RR).*(K1(F).*K2(F).*RR);
H     = 0.5.*(K1(F).*RR + sqrt(Discr))./(1 - RR);
%       if (H <= 0)
%           pHctemp = missingn;
%       else
pHctemp = log(H)./log(0.1);
%       end
varargout{1}=pHctemp;
end % end nested function

function varargout=CalculateTCfrompHfCO2(pHi, fCO2i)
global K0 K1 K2 F;
% ' SUB CalculateTCfrompHfCO2, version 01.02, 12-13-96, written by Ernie Lewis.
% ' Inputs: pH, fCO2, K0, K1, K2
% ' Output: TC
% ' This calculates TC from pH and fCO2, using K0, K1, and K2.
H       = 10.^(-pHi);
TCctemp = K0(F).*fCO2i.*(H.*H + K1(F).*H + K1(F).*K2(F))./(H.*H);
varargout{1}=TCctemp;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCO3 input functions, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculateTAfrompHHCO3(pHi, HCO3i)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculateTAfrompHCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: pH, HCO3, K(), T()
% ' Output: TA
% ' This calculates TA from pH and HCO3.
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
H         = 10.^(-pHi);
CAlk      = HCO3i.*(2.*K2F./H + 1);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHi); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
TActemp     = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
varargout{1}=TActemp;
end % end nested function

function varargout=CalculatepHfromTAHCO3(TAi, HCO3i)
global K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculatepHfromTAHCO3, version 01.0, 8-18, added by J. Sharp with
% ' modifications from Denis Pierrot.
% ' Inputs: TA, CO3, K0, K(), T()
% ' Output: pH
%
% ' This calculates pH from TA and CO3 using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013) and extended by Humphreys et al. (2021).
%
% ' This will continue iterating until all values in the vector are
% ' "abs(deltapH) < pHTol"
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input.
K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl         = sum(F); % vectorlength
% Find initital pH guess using method of Munhoven (2013)
pHGuess    = CalculatepHfromTAHCO3Munhoven(TAi, HCO3i);
ln10       = log(10);
pH         = pHGuess;
pHTol      = 0.0001; % tolerance
deltapH    = pHTol+pH;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    CAlk      = HCO3i.*(H+2.*K2F)./H;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % '               find Slope dTA/dpH
    % '               (this is not exact, but keeps all important terms):
    Slope = ln10 .* (2 .* HCO3i .* K2F ./ H + BAlk .* H ./ (KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculatepHfromTCHCO3(TCi, HCO3i)
global K1 K2 F;
% ' SUB CalculatepHfromTCHCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, HCO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from TC and HCO3 using K1 and K2 by solving the
% '       quadratic in H: TC = HCO3i.*(H./K1 + 1 + K2./H).
% '       Therefore:      0  = H.*H./K1 + (1-TC/HCO3i).*H + K2.
% ' if there is not a real root, then pH is returned as missingn.
RR = TCi./HCO3i;
%       if RR >= 1
%          varargout{1}= missingn;
%          disp('nein!');return;
%       end
% check after sub to see if pH = missingn.
Discr = ((1-RR).*(1-RR) - 4.*(1./(K1(F))).*(K2(F)));
H     = 0.5.*((-(1-RR)) - sqrt(Discr))./(1./(K1(F))); % Subtraction
%       if (H <= 0)
%           pHctemp = missingn;
%       else
pHctemp = log(H)./log(0.1);
%       end
varargout{1}=pHctemp;
end % end nested function

function varargout=CalculatepHfromfCO2HCO3(fCO2i, HCO3i)
global K0 K1 F;
% ' SUB CalculatepHfromfCO2HCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: fCO2, HCO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from fCO2 and HCO3, using K0, K1, and K2.
H            = (fCO2i.*K0(F).*K1(F))./HCO3i;  % removed incorrect (F) index from HCO3i // MPH
pHx          = -log10(H);
varargout{1} = pHx;
end % end nested function

function varargout=CalculatepHfCO2fromTCHCO3(TCx, HCO3x)
% Outputs pH fCO2, in that order
% SUB CalculatepHfCO2fromTCHCO3, version 01.0, 3-19, added by J. Sharp
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, TC, HCO3, Sal, K(), T(), TempC, Pdbar
% Outputs: pH, fCO2
% This calculates pH and fCO2 from TC and HCO3 at output conditions.
pHx   = CalculatepHfromTCHCO3(TCx, HCO3x); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
fCO2x = CalculatefCO2fromTCpH(TCx, pHx);
varargout{1} = pHx;
varargout{2} = fCO2x;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO3 input functions, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculateTAfrompHCO3(pHi, CO3i)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculateTAfrompHCO3, version 01.0, 8-18, added by J. Sharp
% ' Inputs: pH, CO3, K(), T()
% ' Output: TA
% ' This calculates TA from pH and CO3.
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
H         = 10.^(-pHi);
CAlk      = CO3i.*(H./K2F + 2);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHi); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
TActemp     = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
varargout{1}=TActemp;
end % end nested function

function varargout=CalculatepHfromTACO3(TAi, CO3i)
global K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculatepHfromTACO3, version 01.0, 8-18, added by J. Sharp with
% ' modifications from Denis Pierrot.
% ' Inputs: TA, CO3, K0, K(), T()
% ' Output: pH
%
% ' This calculates pH from TA and CO3 using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013) and extended by Humphreys et al. (2021).
%
% ' This will continue iterating until all values in the vector are
% ' "abs(deltapH) < pHTol"
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input.
K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl         = sum(F); % vectorlength
% Find initital pH guess using method of Munhoven (2013)
pHGuess    = CalculatepHfromTACO3Munhoven(TAi, CO3i);
ln10       = log(10);
pH         = pHGuess;
pHTol      = 0.0001; % tolerance
deltapH    = pHTol+pH;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    CAlk      = CO3i.*(H+2.*K2F)./K2F;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % '               find Slope dTA/dpH
    % '               (this is not exact, but keeps all important terms):
    Slope = ln10 .* (-CO3i .* H ./ K2F + BAlk .* H ./ (KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000 
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculatepHfromTCCO3(TCi, CO3i)
global K1 K2 F;
% ' SUB CalculatepHfromTCCO3, version 01.0, 8-18, added by J. Sharp
% ' Inputs: TC, CO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from TC and CO3 using K1 and K2 by solving the
% '       quadratic in H: TC = CO3i.*(H.*H/(K1.*K2) + H./K2 + 1).
% '       Therefore:      0  = H.*H/(K1.*K2) + H./K2 + (1-TC./CO3i).
% ' if there is not a real root, then pH is returned as missingn.
RR = TCi./CO3i;
%       if RR >= 1
%          varargout{1}= missingn;
%          disp('nein!');return;
%       end
% check after sub to see if pH = missingn.
Discr = ((1./K2(F)).*(1./K2(F)) - 4.*(1./(K1(F).*K2(F))).*(1-RR));
H     = 0.5.*((-1./K2(F)) + sqrt(Discr))./(1./(K1(F).*K2(F))); % Addition
%       if (H <= 0)
%           pHctemp = missingn;
%       else
pHctemp = log(H)./log(0.1);
%       end
varargout{1}=pHctemp;
end % end nested function

function varargout=CalculatepHfromfCO2CO3(fCO2i, CO3i)
global K0 K1 K2 F;
% ' SUB CalculatepHfromfCO2CO3, version 01.0, 8-18, added by J. Sharp
% ' Inputs: fCO2, CO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from fCO2 and CO3, using K0, K1, and K2.
H            = sqrt((fCO2i.*K0(F).*K1(F).*K2(F))./CO3i);    % removed incorrect (F) index from CO3i // MPH
pHx          = -log10(H);
varargout{1} = pHx;
end % end nested function

function varargout=CalculatepHfCO2fromTCCO3(TCx, CO3x)
% Outputs pH fCO2, in that order
% SUB CalculatepHfCO2fromTCCO3, version 01.0, 8-18, added by J. Sharp
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, TC, CO3, Sal, K(), T(), TempC, Pdbar
% Outputs: pH, fCO2
% This calculates pH and fCO2 from TC and CO3 at output conditions.
pHx   = CalculatepHfromTCCO3(TCx, CO3x); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
fCO2x = CalculatefCO2fromTCpH(TCx, pHx);
varargout{1} = pHx;
varargout{2} = fCO2x;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO3/HCO3 input function, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculatepHfromCO3HCO3(CO3x, HCO3x)
global K2 F
% ' SUB CalculatepHfromCO3HCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: CO3, HCO3, K2
% ' Output: pH
% ' This calculates fCO2 from TC and pH, using K2.
H            = HCO3x.*K2(F)./CO3x;
pHx          = -log10(H);
varargout{1} = pHx;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO3/HCO3/CO2 output functions, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculateCO3HCO3fromTCpH(TCx, pHx)
global K1 K2 F
% ' SUB CalculateCO3HCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, pH, K1, K2
% ' Output: CO3, HCO3, CO2
% ' This calculates CO3, HCO3, and CO2 from TC and pH, using K1, and K2.
H            = 10.^(-pHx);
CO3x         = TCx.*K1(F).*K2(F)./(K1(F).*H + H.*H + K1(F).*K2(F));
HCO3x        = TCx.*K1(F).*H./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = CO3x;
varargout{2} = HCO3x;
end % end nested function

function varargout=CalculateCO3fromTCpH(TCx, pHx)
global K1 K2 F
% ' SUB CalculateCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, pH, K1, K2
% ' Output: CO3, CO2
% ' This calculates CO3 and CO2 from TC and pH, using K1, and K2.
H            = 10.^(-pHx);
CO3x         = TCx.*K1(F).*K2(F)./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = CO3x;
end % end nested function

function varargout=CalculateHCO3fromTCpH(TCx, pHx)
global K1 K2 F
% ' SUB CalculateHCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, pH, K1, K2
% ' Output: HCO3, CO2
% ' This calculates HCO3 and CO2 from TC and pH, using K1, and K2.
H            = 10.^(-pHx);
HCO3x        = TCx.*K1(F).*H./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = HCO3x;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial pH estimates via Munhoven (2013), Humphreys et al (2021).
% Added by J. Sharp (3-30-21)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculatepHfromTATCMunhoven(TAi, TCi)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = K1F.*K2F.*KBF.*(1-(2.*TCi+TBF)./TAi);
g1 = K1F.*(KBF.*(1-TBF./TAi-TCi./TAi)+K2F.*(1-2.*TCi./TAi));
g2 = KBF.*(1-TBF./TAi)+K1F.*(1-TCi./TAi);
% Determine g21min
g21min = g2.^2-3.*g1;
g21min_positive = g21min > 0;
sq21 = nan(size(TAi,1),1);
sq21(g21min_positive) = sqrt(g21min(g21min_positive));
sq21(~g21min_positive) = 0;
% Determine Hmin
Hmin = nan(size(TAi,1),1);
g2_positive = g2 >=0;
Hmin(~g2_positive) = (-g2(~g2_positive) + sq21(~g2_positive))./3;
Hmin(g2_positive) = -g1(g2_positive)./(g2(g2_positive) + sq21(g2_positive));
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= 0;
pHGuess(idx) = -log10(1e-3);
idx = TAi > 0 & TAi < 2.*TCi + TBF;
pHGuess(idx & g21min_positive) = ...
    -log10(Hmin(idx & g21min_positive) + ...
    sqrt(-(Hmin(idx & g21min_positive).^3 + g2(idx & g21min_positive).*Hmin(idx & g21min_positive).^2 + ...
    g1(idx & g21min_positive).*Hmin(idx & g21min_positive) + ...
    g0(idx & g21min_positive))./sq21(idx & g21min_positive)));
pHGuess(idx & ~g21min_positive) = -log10(1e-7);
idx = TAi >= 2.*TCi + TBF;
pHGuess(idx) = -log10(1e-10);
varargout{1}=pHGuess;
end

function varargout=CalculatepHfromTACO2Munhoven(TAi, CO2x)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = -2.*K1F.*K2F.*KBF.*CO2x./TAi;
g1 = -K1F.*(2.*K2F.*CO2x+KBF.*CO2x)./TAi;
g2 = KBF-(TBF.*KBF+K1F.*CO2x)./TAi;
% Determine Hmin
g21min = g2.^2-3.*g1;
g21min_positive = g21min > 0;
sq21 = nan(size(TAi,1),1);
sq21(g21min_positive) = sqrt(g21min(g21min_positive));
sq21(~g21min_positive) = 0;
Hmin = nan(size(TAi,1),1);
g2_positive = g2 >=0;
Hmin(~g2_positive) = (-g2(~g2_positive) + sq21(~g2_positive))./3;
Hmin(g2_positive) = -g1(g2_positive)./(g2(g2_positive) + sq21(g2_positive));
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= 0;
pHGuess(idx) = -log10(1e-3);
idx = TAi > 0;
pHGuess(idx & g21min_positive) = ...
    -log10(Hmin(idx & g21min_positive) + ...
    sqrt(-(Hmin(idx & g21min_positive).^3 + g2(idx & g21min_positive).*Hmin(idx & g21min_positive).^2 + ...
    g1(idx & g21min_positive).*Hmin(idx & g21min_positive)+...
    g0(idx & g21min_positive))./sq21(idx & g21min_positive)));
pHGuess(idx & ~g21min_positive) = -log10(1e-7);
varargout{1}=pHGuess;
end

function varargout=CalculatepHfromTAHCO3Munhoven(TAi, HCO3x)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = 2.*K2F.*KBF.*HCO3x;
g1 = KBF.*(HCO3x+TBF-TAi)+2.*K2F.*HCO3x;
g2 = HCO3x-TAi;
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= HCO3x;
pHGuess(idx) = -log10(1e-3);
idx = TAi > HCO3x;
pHGuess(idx) = ...
    -log10((-g1(idx)-sqrt(g1(idx).^2-4.*g0(idx).*g2(idx)))./(2.*g2(idx)));
varargout{1}=pHGuess;
end

function varargout=CalculatepHfromTACO3Munhoven(TAi, CO3x)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = K2F.*KBF.*(2.*CO3x+TBF-TAi);
g1 = KBF.*CO3x+K2F.*(2.*CO3x-TAi);
g2 = CO3x;
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= 2.*CO3x+TBF;
pHGuess(idx) = -log10(1e-3);
idx = TAi > 2.*CO3x+TBF;
pHGuess(idx) = ...
    -log10((-g1(idx)+sqrt(g1(idx).^2-4.*g0(idx).*g2(idx)))./(2.*g2(idx)));
varargout{1}=pHGuess;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=RevelleFactor(TAi, TCi)
% global WhichKs;
% ' SUB RevelleFactor, version 01.03, 01-07-97, written by Ernie Lewis.
% ' Inputs: WhichKs%, TA, TC, K0, K(), T()
% ' Outputs: Revelle
% ' This calculates the Revelle factor (dfCO2/dTC)|TA/(fCO2/TC).
% ' It only makes sense to talk about it at pTot = 1 atm, but it is computed
% '       here at the given K(), which may be at pressure <> 1 atm. Care must
% '       thus be used to see if there is any validity to the number computed.
TC0 = TCi;
dTC = 0.00000001;% ' 0.01 umol/kg-SW (lower than prior versions of CO2SYS)
% ' Find fCO2 at TA, TC + dTC
TCi = TC0 + dTC;
pHc= CalculatepHfromTATC(TAi, TCi);
fCO2c= CalculatefCO2fromTCpH(TCi, pHc);
fCO2plus = fCO2c;
% ' Find fCO2 at TA, TC - dTC
TCi = TC0 - dTC;
pHc= CalculatepHfromTATC(TAi, TCi);
fCO2c= CalculatefCO2fromTCpH(TCi, pHc);
fCO2minus = fCO2c;
% CalculateRevelleFactor:
Revelle = (fCO2plus - fCO2minus)./dTC./((fCO2plus + fCO2minus)./TC0); % Corrected error pointed out by MP Humphreys (https://pyco2sys.readthedocs.io/en/latest/validate/)
varargout{1}=Revelle;
end % end nested function


function varargout=CalculateAlkParts(pH)
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F;
% ' SUB CalculateAlkParts, version 01.03, 10-10-97, written by Ernie Lewis.
% ' Inputs: pH, TC, K(), T()
% ' Outputs: BAlk, OH, PAlk, SiAlk, Hfree, HSO4, HF
% ' This calculates the various contributions to the alkalinity.
% ' Though it is coded for H on the total pH scale, for the pH values occuring
% ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
% ' negligible) as long as the K Constants are on that scale.

KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);

H         = 10.^(-pH);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree); %' since KF is on the free scale

varargout{1} = BAlk;  varargout{2} = OH; varargout{3} = PAlk;
varargout{4} = SiAlk; varargout{5} = AmmAlk; varargout{6} = HSAlk;
varargout{7} = Hfree; varargout{8} = HSO4; varargout{9} = HF;
end % end nested function


function varargout=CaSolubility(Sal, TempC, Pdbar, TC, pH)
global K1 K2 TempK logTempK sqrSal Pbar RT WhichKs CAL ntps F
global PertK    % Id of perturbed K
global Perturb  % perturbation
% '***********************************************************************
% ' SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
% ' Inputs: WhichKs%, Sal, TempCi, Pdbari, TCi, pHi, K1, K2
% ' Outputs: OmegaCa, OmegaAr
% ' This calculates omega, the solubility ratio, for calcite and aragonite.
% ' This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
% '       where Ksp is the solubility product (either KCa or KAr).
% '***********************************************************************
% ' These are from:
% ' Mucci, Alphonso, The solubility of calcite and aragonite in seawater
% '       at various salinities, temperatures, and one atmosphere total
% '       pressure, American Journal of Science 283:781-799, 1983.
% ' Ingle, S. E., Solubility of calcite in the ocean,
% '       Marine Chemistry 3:301-319, 1975,
% ' Millero, Frank, The thermodynamics of the carbonate system in seawater,
% '       Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
% ' Ingle et al, The solubility of calcite in seawater at atmospheric pressure
% '       and 35%o salinity, Marine Chemistry 1:295-307, 1973.
% ' Berner, R. A., The solubility of calcite and aragonite in seawater in
% '       atmospheric pressure and 34.5%o salinity, American Journal of
% '       Science 276:713-730, 1976.
% ' Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
% ' Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
% '       boric acid, and the pHi of seawater, Limnology and Oceanography
% '       13:403-417, 1968.
% '***********************************************************************
Ca=CAL(F);
Ar=nan(sum(F),1);
KCa=nan(sum(F),1);
KAr=nan(sum(F),1);
TempKx=TempK(F);
logTempKx=logTempK(F);
sqrSalx=sqrSal(F);
Pbarx=Pbar(F);
RTx=RT(F);
FF=(WhichKs(F)~=6 & WhichKs(F)~=7);
if any(FF)
% (below here, F isn't used, since almost always all rows match the above criterium,
%  in all other cases the rows will be overwritten later on).
    % CalciteSolubility:
    % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKCa = -171.9065 - 0.077993.*TempKx(FF) + 2839.319./TempKx(FF);
    logKCa = logKCa + 71.595.*logTempKx(FF)./log(10);
    logKCa = logKCa + (-0.77712 + 0.0028426.*TempKx(FF) + 178.34./TempKx(FF)).*sqrSalx(FF);
    logKCa = logKCa - 0.07711.*Sal(FF) + 0.0041249.*sqrSalx(FF).*Sal(FF);
    % '       sd fit = .01 (for Sal part, not part independent of Sal)
    KCa(FF) = 10.^(logKCa);% ' this is in (mol/kg-SW)^2
    % AragoniteSolubility:
    % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKAr = -171.945 - 0.077993.*TempKx(FF) + 2903.293./TempKx(FF);
    logKAr = logKAr + 71.595.*logTempKx(FF)./log(10);
    logKAr = logKAr + (-0.068393 + 0.0017276.*TempKx(FF) + 88.135./TempKx(FF)).*sqrSalx(FF);
    logKAr = logKAr - 0.10018.*Sal(FF) + 0.0059415.*sqrSalx(FF).*Sal(FF);
    % '       sd fit = .009 (for Sal part, not part independent of Sal)
    KAr(FF)    = 10.^(logKAr);% ' this is in (mol/kg-SW)^2
    % PressureCorrectionForCalcite:
    % '       Ingle, Marine Chemistry 3:301-319, 1975
    % '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
    % '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
    deltaVKCa = -48.76 + 0.5304.*TempC(FF);
    KappaKCa  = (-11.76 + 0.3692.*TempC(FF))./1000;
    lnKCafac  = (-deltaVKCa + 0.5.*KappaKCa.*Pbarx(FF)).*Pbarx(FF)./RTx(FF);
    KCa(FF)       = KCa(FF).*exp(lnKCafac);
    % PressureCorrectionForAragonite:
    % '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
    % '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
    % '       and 10^3 for Kappa factor)
    deltaVKAr = deltaVKCa + 2.8;
    KappaKAr  = KappaKCa;
    lnKArfac  = (-deltaVKAr + 0.5.*KappaKAr.*Pbarx(FF)).*Pbarx(FF)./RTx(FF);
    KAr(FF)       = KAr(FF).*exp(lnKArfac);
end
FF=(WhichKs(F)==6 | WhichKs(F)==7);
if any(FF)
    % *** CalculateKCaforGEOSECS:
    % Ingle et al, Marine Chemistry 1:295-307, 1973 is referenced in
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
    % but the fit is actually from Ingle, Marine Chemistry 3:301-319, 1975)
    KCa(FF) = 0.0000001.*(-34.452 - 39.866.*Sal(FF).^(1./3) +...
        110.21.*log(Sal(FF))./log(10) - 0.0000075752.*TempKx(FF).^2);
    % this is in (mol/kg-SW)^2
    %
    % *** CalculateKArforGEOSECS:
    % Berner, R. A., American Journal of Science 276:713-730, 1976:
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
    KAr(FF) = 1.45.*KCa(FF);% ' this is in (mol/kg-SW)^2
    % Berner (p. 722) states that he uses 1.48.
    % It appears that 1.45 was used in the GEOSECS calculations
    %
    % *** CalculatePressureEffectsOnKCaKArGEOSECS:
    % Culberson and Pytkowicz, Limnology and Oceanography 13:403-417, 1968
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
    % but their paper is not even on this topic).
    % The fits appears to be new in the GEOSECS report.
    % I can't find them anywhere else.
    KCa(FF) = KCa(FF).*exp((36   - 0.2 .*TempC(FF)).*Pbarx(FF)./RTx(FF));
    KAr(FF) = KAr(FF).*exp((33.3 - 0.22.*TempC(FF)).*Pbarx(FF)./RTx(FF));
end
% Added by JM Epitalon
% For computing derivative with respect to KCa or KAr, one has to perturb the value of one K
% Requested perturbation is passed through global variables PertK and Perturb
if (~ isempty(PertK))
    switch PertK
        case {'KSPA'}   % solubility Product for Aragonite
            KAr = KAr + Perturb;
        case {'KSPC'}   % for Calcite
            KCa = KCa + Perturb;
        case {'CAL'}   % for calcium concentration
            Ca  = Ca  + Perturb;
    end
end

% CalculateOmegasHere:
H = 10.^(-pH);
CO3 = TC.*K1(F).*K2(F)./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = CO3.*Ca./KCa; % OmegaCa, dimensionless
varargout{2} = CO3.*Ca./KAr; % OmegaAr, dimensionless
end % end nested function

function varargout=FindpHOnAllScales(pH)
global pHScale K T TS KS TF KF fH F ntps;
% ' SUB FindpHOnAllScales, version 01.02, 01-08-97, written by Ernie Lewis.
% ' Inputs: pHScale%, pH, K(), T(), fH
% ' Outputs: pHNBS, pHfree, pHTot, pHSWS
% ' This takes the pH on the given scale and finds the pH on all scales.
%  TS = T(3); TF = T(2);
%  KS = K(6); KF = K(5);% 'these are at the given T, S, P
TSx=TS(F); KSx=KS(F); TFx=TF(F); KFx=KF(F);fHx=fH(F);
FREEtoTOT = (1 + TSx./KSx); % ' pH scale conversion factor
SWStoTOT  = (1 + TSx./KSx)./(1 + TSx./KSx + TFx./KFx);% ' pH scale conversion factor
factor=nan(sum(F),1);
nF=pHScale(F)==1;  %'"pHtot"
factor(nF) = 0;
nF=pHScale(F)==2; % '"pHsws"
factor(nF) = -log(SWStoTOT(nF))./log(0.1);
nF=pHScale(F)==3; % '"pHfree"
factor(nF) = -log(FREEtoTOT(nF))./log(0.1);
nF=pHScale(F)==4;  %'"pHNBS"
factor(nF) = -log(SWStoTOT(nF))./log(0.1) + log(fHx(nF))./log(0.1);
pHtot  = pH    - factor;    % ' pH comes into this sub on the given scale
pHNBS  = pHtot - log(SWStoTOT) ./log(0.1) + log(fHx)./log(0.1);
pHfree = pHtot - log(FREEtoTOT)./log(0.1);
pHsws  = pHtot - log(SWStoTOT) ./log(0.1);
varargout{1}=pHtot;
varargout{2}=pHsws;
varargout{3}=pHfree;
varargout{4}=pHNBS;
end % end nested function
