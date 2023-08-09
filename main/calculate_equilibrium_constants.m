
function calculate_equilibrium_constants(TempC,Pdbar)
global pHScale WhichKs WhoseKSO4 WhoseKF WhoseTB sqrSal Pbar RT;
global K0 fH FugFac VPFac ntps TempK logTempK;
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi CAL RGasConstant Sal p_opt;

% SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar
% Outputs: K0, K(), T(), fH, FugFac, VPFac
% This finds the Constants of the CO2 system in seawater or freshwater,
% corrects them for pressure, and reports them on the chosen pH scale.
% The process is as follows: the Constants (except KS, KF which stay on the
% free scale - these are only corrected for pressure) are
%       1) evaluated as they are given in the literature
%       2) converted to the SWS scale in mol/kg-SW or to the NBS scale
%       3) corrected for pressure
%       4) converted to the SWS pH scale in mol/kg-SW
%       5) converted to the chosen pH scale
%
%       PROGRAMMER'S NOTE: all logs are log base e
%       PROGRAMMER'S NOTE: all Constants are converted to the pH scale
%               pHScale% (the chosen one) in units of mol/kg-SW
%               except KS and KF are on the free scale
%               and KW is in units of (mol/kg-SW)^2
TempK    = TempC + 273.15;
RT       = RGasConstant.*TempK;
logTempK = log(TempK);
Pbar     = Pdbar ./ 10;

% Generate empty vectors for holding results
TB = nan(ntps,1);
TF = nan(ntps,1);
TS = nan(ntps,1);
CAL = nan(ntps,1);

% CalculateTB - Total Borate:
F=(WhichKs==8); % Pure water case.
if any(F)
    TB(F) = 0;
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    TB(F) = 0.0004106.*Sal(F)./35; % in mol/kg-SW
    % this is .00001173.*Sali
    % this is about 1% lower than Uppstrom's value
    % Culkin, F., in Chemical Oceanography,
    % ed. Riley and Skirrow, 1965:
    % GEOSECS references this, but this value is not explicitly
    % given here
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8); % All other cases
if any(F)
	FF=F&(WhoseTB==1); % If user opted for Uppstrom's values:
	if any(FF)
	    % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
	    % this is .000416.*Sali./35. = .0000119.*Sali
		% TB(FF) = (0.000232./10.811).*(Sal(FF)./1.80655); % in mol/kg-SW
	    TB(FF) =  0.0004157.*Sal(FF)./35; % in mol/kg-SW
	end
	FF=F&(WhoseTB==2); % If user opted for the Lee et al. values:
	if any(FF)
		% Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.	
	 	% Geochimica Et Cosmochimica Acta 74 (6): 1801-1811.
		TB(FF) =  0.0004326.*Sal(FF)./35; % in mol/kg-SW
	end
end

% CalculateCAL - Total Calcium:
F=(WhichKs~=6 & WhichKs~=7);
    % Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    % this is .010285.*Sali./35
    CAL(F) = 0.02128./40.087.*(Sal(F)./1.80655);
F=(WhichKs==6 | WhichKs==7);
    % *** CalculateCaforGEOSECS:
    % Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982) 
    CAL(F) = 0.01026.*Sal(F)./35;

% CalculateTF;
% Riley, J. P., Deep-Sea Research 12:219-220, 1965:
% this is .000068.*Sali./35. = .00000195.*Sali
TF = (0.000067./18.998).*(Sal./1.80655); % in mol/kg-SW

% CalculateTS ;
% Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
% this is .02824.*Sali./35. = .0008067.*Sali
TS = (0.14./96.062).*(Sal./1.80655); % in mol/kg-SW

% CalculateK0:
% Weiss, R. F., Marine Chemistry 2:203-215, 1974.
TempK100  = TempK./100;
lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + Sal .*...
    (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
K0   = exp(lnK0);                 % this is in mol/kg-SW/atm
vCO2 = 32.3;                      % partial molal volume of CO2 (cm3 / mol)
                                  % from Weiss (1974, Appendix, paragraph 3)
if p_opt == 0
    prscorr = 1; % Set pressure correction to 1
elseif p_opt == 1
    prscorr = exp((-Pbar).*vCO2./RT); % Calculate pressure correction to K0
else
    disp('co2_press must be set to either 0 or 1'); % Display error message
end         
K0   = K0 .* prscorr; % this is in mol/kg-SW/atm

% CalculateIonS:
% This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
IonS         = 19.924 .* Sal ./ (1000 - 1.005   .* Sal);

% CalculateKS:
lnKS   = nan(ntps,1); pKS  = nan(ntps,1); KS   = nan(ntps,1);
logKS0 = nan(ntps,1); logKSK0 = nan(ntps,1);
F=(WhoseKSO4==1);
if any(F)
    % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
    % The goodness of fit is .021.
    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    % TYPO on p. 121: the constant e9 should be e8.
    % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
  lnKS(F) = -4276.1./TempK(F) + 141.328 - 23.093.*logTempK(F) +...             
      (-13856./TempK(F) + 324.57 - 47.986.*logTempK(F)).*sqrt(IonS(F)) +...     
      (35474./TempK(F) - 771.54 + 114.723.*logTempK(F)).*IonS(F) +...           
      (-2698./TempK(F)).*sqrt(IonS(F)).*IonS(F) + (1776./TempK(F)).*IonS(F).^2; 
	KS(F) = exp(lnKS(F))...            % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005 .* Sal(F));   % convert to mol/kg-SW
end
F=(WhoseKSO4==2);
if any(F)
    % Khoo et al, Analytical Chemistry, 49(1):29-34, 1977
    % KS was found by titrations with a hydrogen electrode
    % of artificial seawater containing sulfate (but without F)
    % at 3 salinities from 20 to 45 and artificial seawater NOT
    % containing sulfate (nor F) at 16 salinities from 15 to 45,
    % both at temperatures from 5 to 40 deg C.
    % KS is on the Free pH scale (inherently so).
    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    % He finds log(beta) which = my pKS;
    % his beta is an association constant.
    % The rms error is .0021 in pKS, or about .5% in KS.
    % This is equation 20 on p. 33:
    pKS(F) = 647.59 ./ TempK(F) - 6.3451 + 0.019085.*TempK(F) - 0.5208.*sqrt(IonS(F));
    KS(F) = 10.^(-pKS(F))...          % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005.*Sal(F));    % convert to mol/kg-SW
end
F=(WhoseKSO4==3);
if any(F)
    % Waters and Millero, Marine Chemistry, 149: 8-22, 2013, with corrections from
    % Waters et al, Marine Chemistry, 165: 66-67, 2014
    logKS0(F) = 562.69486 - 102.5154.*logTempK(F) - 0.0001117033.*TempK(F).*TempK(F) + ...
        0.2477538.*TempK(F) - 13273.76./TempK(F);
    logKSK0(F) = (4.24666 - 0.152671.*TempK(F) + 0.0267059.*TempK(F).*logTempK(F) - 0.000042128.*TempK(F).*TempK(F)).*Sal(F).^0.5 + ...
        (0.2542181 - 0.00509534.*TempK(F) + 0.00071589.*TempK(F).*logTempK(F)).*Sal(F) + (-0.00291179 + 0.0000209968.*TempK(F)).*Sal(F).^1.5 + ...
        -0.0000403724.*Sal(F).^2;
    KS(F) = ((10.^(logKSK0(F))).*(10.^logKS0(F))) ... % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005.*Sal(F));                    % convert to mol/kg-SW
end

% CalculateKF:
KF = NaN(ntps, 1);  % added preallocation here and F-indexing below // MPH
F=(WhoseKF==1);
if any(F)
    % Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
    lnKF = 1590.2./TempK - 12.641 + 1.525.*IonS.^0.5;
    KF(F)   = exp(lnKF(F))...                 % this is on the free pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F));          % convert to mol/kg-SW
end
F=(WhoseKF==2);
if any(F)
    % Perez and Fraga 1987 (to be used for S: 10-40, T: 9-33)
    % P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
    lnKF = 874./TempK - 9.68 + 0.111.*Sal.^0.5;
    KF(F)   = exp(lnKF(F));                   % this is on the free pH scale in mol/kg-SW
end

% CalculatepHScaleConversionFactors:
%       These are NOT pressure-corrected.
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

% CalculatefH
fH = nan(ntps,1);
% Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
F=(WhichKs==8);
if any(F)
    fH(F) = 1; % this shouldn't occur in the program for this case
end
F=(WhichKs==7);
if any(F)
    fH(F) = 1.29 - 0.00204.*  TempK(F) + (0.00046 -...
        0.00000148.*TempK(F)).*Sal(F).*Sal(F);
    % Peng et al, Tellus 39B:439-458, 1987:
    % They reference the GEOSECS report, but round the value
    % given there off so that it is about .008 (1%) lower. It
    % doesn't agree with the check value they give on p. 456.
end
F=(WhichKs~=7 & WhichKs~=8);
if any(F)
    fH(F) = 1.2948 - 0.002036.*TempK(F) + (0.0004607 -...
        0.000001475.*TempK(F)).*Sal(F).^2;
    % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
    % v. 3, 1982 (p. 80);
end

% CalculateKB:
KB      = nan(ntps,1); logKB   = nan(ntps,1);
lnKBtop = nan(ntps,1); lnKB    = nan(ntps,1);
F=(WhichKs==8); % Pure water case
if any(F)
    KB(F) = 0;
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    % This is for GEOSECS and Peng et al.
    % Lyman, John, UCLA Thesis, 1957
    % fit by Li et al, JGR 74:5507-5525, 1969:
    logKB(F) = -9.26 + 0.00886.*Sal(F) + 0.01.*TempC(F);
    KB(F) = 10.^(logKB(F))...  % this is on the NBS scale
        ./fH(F);               % convert to the SWS scale
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
    lnKBtop(F) = -8966.9 - 2890.53.*sqrSal(F) - 77.942.*Sal(F) +...
        1.728.*sqrSal(F).*Sal(F) - 0.0996.*Sal(F).^2;
    lnKB(F) = lnKBtop(F)./TempK(F) + 148.0248 + 137.1942.*sqrSal(F) +...
        1.62142.*Sal(F) + (-24.4344 - 25.085.*sqrSal(F) - 0.2474.*...
        Sal(F)).*logTempK(F) + 0.053105.*sqrSal(F).*TempK(F);
    KB(F) = exp(lnKB(F))...    % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);         % convert to SWS pH scale
end

% CalculateKW:
lnKW = nan(ntps,1); KW = nan(ntps,1);
F=(WhichKs==7);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F) +...
        (-79.2447 + 3298.72./TempK(F) + 12.0408.*logTempK(F)).*...
        sqrSal(F) - 0.019813.*Sal(F);
end
F=(WhichKs==8);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    % refit data of Harned and Owen, The Physical Chemistry of
    % Electrolyte Solutions, 1958
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F);
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
    % his check value of 1.6 umol/kg-SW should be 6.2
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F) +...
        (-5.977 + 118.67./TempK(F) + 1.0495.*logTempK(F)).*...
        sqrSal(F) - 0.01615.*Sal(F);
end
KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
F=(WhichKs==6);
if any(F)
    KW(F) = 0; % GEOSECS doesn't include OH effects
end

% CalculateKP1KP2KP3KSi:
KP1      = nan(ntps,1); KP2      = nan(ntps,1);
KP3      = nan(ntps,1); KSi      = nan(ntps,1);
lnKP1    = nan(ntps,1); lnKP2    = nan(ntps,1);
lnKP3    = nan(ntps,1); lnKSi    = nan(ntps,1);
F=(WhichKs==7);
if any(F)
    KP1(F) = 0.02;
    % Peng et al don't include the contribution from this term,
    % but it is so small it doesn't contribute. It needs to be
    % kept so that the routines work ok.
    % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
    % Limnology and Oceanography 12:243-252, 1967:
    % these are only for sals 33 to 36 and are on the NBS scale
    KP2(F) = exp(-9.039 - 1450./TempK(F))... % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    KP3(F) = exp(4.466 - 7276./TempK(F))...  % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    % Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
    % The Chemical Society (London), Special Publ. 17:751, 1964:
    KSi(F) = 0.0000000004...              % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    
end
F=(WhichKs==6 | WhichKs==8);
if any(F)
    KP1(F) = 0; KP2(F) = 0; KP3(F) = 0; KSi(F) = 0;
    % Neither the GEOSECS choice nor the freshwater choice
    % include contributions from phosphate or silicate.
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
    % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
    % KSi was given on the SWS pH scale in molal units.
    lnKP1(F) = -4576.752./TempK(F) + 115.54 - 18.453.*logTempK(F) + (-106.736./TempK(F) +...
        0.69171).*sqrSal(F) + (-0.65643./TempK(F) - 0.01844).*Sal(F);
    KP1(F) = exp(lnKP1(F));
    lnKP2(F) = -8814.715./TempK(F) + 172.1033 - 27.927.*logTempK(F) + (-160.34./TempK(F) +...
        1.3566).*sqrSal(F) + (0.37335./TempK(F) - 0.05778).*Sal(F);
    KP2(F) = exp(lnKP2(F));
    lnKP3(F) = -3070.75./TempK(F) - 18.126 + (17.27039./TempK(F) + 2.81197).*sqrSal(F) +...
        (-44.99486./TempK(F) - 0.09984).*Sal(F);
    KP3(F) = exp(lnKP3(F));
    lnKSi(F) = -8904.2./TempK(F) + 117.4 - 19.334.*logTempK(F) + (-458.79./TempK(F) +...
        3.5913).*sqrt(IonS(F)) + (188.74./TempK(F) - 1.5998).*IonS(F) +...
        (-12.1652./TempK(F) + 0.07871).*IonS(F).^2;
    KSi(F) = exp(lnKSi(F))...                % this is on the SWS pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F));        % convert to mol/kg-SW
end

% Calculate KNH4 and KH2S: added by J. Sharp
KNH4           = nan(ntps,1); KH2S       = nan(ntps,1);
PKNH4expCW     = nan(ntps,1); lnKH2S     = nan(ntps,1);
F=(WhichKs==6 | WhichKs==7 | WhichKs==8); % GEOSECS or freshwater cases
if any(F)
    KNH4(F) = 0;
    KH2S(F) = 0;
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8); % All other cases
if any(F)
% Ammonia dissociation constant from Yao and Millero (1995)
%   KNH4(F) = (exp(-6285.33./TempK(F)+0.0001635.*TempK(F)-0.25444+...
%             (0.46532-123.7184./TempK(F)).*Sal(F).^0.5+(-0.01992+...
%             3.17556./TempK(F)).*Sal(F)))...
%             ./SWStoTOT(F);                    % convert to SWS pH scale
% Ammonia dissociation constant from Clegg and Whitfield (1995)
  PKNH4(F) = 9.244605-2729.33.*(1./298.15-1./TempK(F)) +...
          (0.04203362-11.24742./TempK(F)).*Sal(F).^0.25+...  % added missing (F) index on Sal // MPH
          (-13.6416+1.176949.*TempK(F).^0.5-...
          0.02860785.*TempK(F)+545.4834./TempK(F)).*Sal(F).^0.5+...
          (-0.1462507+0.0090226468.*TempK(F).^0.5-...
          0.0001471361.*TempK(F)+10.5425./TempK(F)).*Sal(F).^1.5+...
          (0.004669309-0.0001691742.*TempK(F).^0.5-...
          0.5677934./TempK(F)).*Sal(F).^2+...
          (-2.354039E-05+0.009698623./TempK(F)).*Sal(F).^2.5;
  KNH4(F)  = 10.^-PKNH4(F);                    % total scale, mol/kg-H2O
  KNH4(F)  = KNH4(F).*(1-0.001005.*Sal(F)); % mol/kg-SW
  KNH4(F)  = KNH4(F)./SWStoTOT(F);             % converts to SWS pH scale

% First hydrogen sulfide dissociation constant from Millero et al. (1988)
  KH2S(F)  = (exp(225.838-13275.3./TempK(F)-34.6435.*log(TempK(F))+...
              0.3449.*Sal(F).^0.5-0.0274.*Sal(F)))...
              ./SWStoTOT(F);                    % convert to SWS pH scale

end

% CalculateK1K2:
logK1    = nan(ntps,1); lnK1     = nan(ntps,1);
pK1      = nan(ntps,1); K1       = nan(ntps,1);
logK2    = nan(ntps,1); lnK2     = nan(ntps,1);
pK2      = nan(ntps,1); K2       = nan(ntps,1);
F=(WhichKs==1);
if any(F)
    % ROY et al, Marine Chemistry, 44:249-267, 1993
    % (see also: Erratum, Marine Chemistry 45:337, 1994
    % and Erratum, Marine Chemistry 52:183, 1996)
    % Typo: in the abstract on p. 249: in the eq. for lnK1* the
    % last term should have S raised to the power 1.5.
    % They claim standard deviations (p. 254) of the fits as
    % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
    % They also claim (p. 258) 2s precisions of .004 in pK1 and
    % .006 in pK2. These are consistent, but Andrew Dickson
    % (personal communication) obtained an rms deviation of about
    % .004 in pK1 and .003 in pK2. This would be a 2s precision
    % of about 2% in K1 and 1.5% in K2.
    % T:  0-45  S:  5-45. Total Scale. Artificial sewater.
    % This is eq. 29 on p. 254 and what they use in their abstract:
    lnK1(F) = 2.83655 - 2307.1266./TempK(F) - 1.5529413.*logTempK(F) +...
        (-0.20760841 - 4.0484./TempK(F)).*sqrSal(F) + 0.08468345.*Sal(F) -...
        0.00654208.*sqrSal(F).*Sal(F);
    K1(F) = exp(lnK1(F))...            % this is on the total pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F))...    % convert to mol/kg-SW
        ./SWStoTOT(F);                 % convert to SWS pH scale
    % This is eq. 30 on p. 254 and what they use in their abstract:
    lnK2(F) = -9.226508 - 3351.6106./TempK(F) - 0.2005743.*logTempK(F) +...
        (-0.106901773 - 23.9722./TempK(F)).*sqrSal(F) + 0.1130822.*Sal(F) -...
        0.00846934.*sqrSal(F).*Sal(F);
    K2(F) = exp(lnK2(F))...            % this is on the total pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F))...    % convert to mol/kg-SW
        ./SWStoTOT(F);                 % convert to SWS pH scale
end
F=(WhichKs==2);
if any(F)
    % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
    % The 2s precision in pK1 is .011, or 2.5% in K1.
    % The 2s precision in pK2 is .02, or 4.5% in K2.
    % This is in Table 5 on p. 1652 and what they use in the abstract:
    pK1(F) = 812.27./TempK(F) + 3.356 - 0.00171.*Sal(F).*logTempK(F)...
        + 0.000091.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 5 on p. 1652 and what they use in the abstract:
    pK2(F) = 1450.87./TempK(F) + 4.604 - 0.00385.*Sal(F).*logTempK(F)...
        + 0.000182.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==3);
if any(F)
    % HANSSON refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
    % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
    % on the SWS pH scale in mol/kg-SW.
    % Hansson gave his results on the Total scale (he called it
    % the seawater scale) and in mol/kg-SW.
    % Typo in DM on p. 1739 in Table 4: the equation for pK2*
    % for Hansson should have a .000132 *S^2
    % instead of a .000116 *S^2.
    % The 2s precision in pK1 is .013, or 3% in K1.
    % The 2s precision in pK2 is .017, or 4.1% in K2.
    % This is from Table 4 on p. 1739.
    pK1(F) = 851.4./TempK(F) + 3.237 - 0.0106.*Sal(F) + 0.000105.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is from Table 4 on p. 1739.
    pK2(F) = -3885.4./TempK(F) + 125.844 - 18.141.*logTempK(F)...
        - 0.0192.*Sal(F) + 0.000132.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==4);
if any(F)
    % MEHRBACH refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
    % on the SWS pH scale in mol/kg-SW.
    % Mehrbach et al gave results on the NBS scale.
    % The 2s precision in pK1 is .011, or 2.6% in K1.
    % The 2s precision in pK2 is .020, or 4.6% in K2.
	% Valid for salinity 20-40.
    % This is in Table 4 on p. 1739.
    pK1(F) = 3670.7./TempK(F) - 62.008 + 9.7944.*logTempK(F)...
             - 0.0118.*Sal(F) + 0.000116.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 4 on p. 1739.
    pK2(F) = 1394.7./TempK(F) + 4.777 - 0.0184.*Sal(F) + 0.000118.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==5);
if any(F)
    % HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
    % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
    % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
    % on the SWS pH scale in mol/kg-SW.
    % Typo in DM on p. 1740 in Table 5: the second equation
    % should be pK2* =, not pK1* =.
    % The 2s precision in pK1 is .017, or 4% in K1.
    % The 2s precision in pK2 is .026, or 6% in K2.
	% Valid for salinity 20-40.
    % This is in Table 5 on p. 1740.
    pK1(F) = 845./TempK(F) + 3.248 - 0.0098.*Sal(F) + 0.000087.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 5 on p. 1740.
    pK2(F) = 1377.3./TempK(F) + 4.824 - 0.0185.*Sal(F) + 0.000122.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
    % Limnology and Oceanography, 18(6):897-907, 1973.
	% I.e., these are the original Mehrbach dissociation constants.
    % The 2s precision in pK1 is .005, or 1.2% in K1.
    % The 2s precision in pK2 is .008, or 2% in K2.
    pK1(F) = - 13.7201 + 0.031334.*TempK(F) + 3235.76./TempK(F)...
        + 1.3e-5*Sal(F).*TempK(F) - 0.1032.*Sal(F).^0.5;
    K1(F) = 10.^(-pK1(F))...         % this is on the NBS scale
        ./fH(F);                     % convert to SWS scale
    pK2(F) = 5371.9645 + 1.671221.*TempK(F) + 0.22913.*Sal(F) + 18.3802.*log10(Sal(F))...
             - 128375.28./TempK(F) - 2194.3055.*log10(TempK(F)) - 8.0944e-4.*Sal(F).*TempK(F)...
             - 5617.11.*log10(Sal(F))./TempK(F) + 2.136.*Sal(F)./TempK(F); % pK2 is not defined for Sal=0, since log10(0)=-inf
    K2(F) = 10.^(-pK2(F))...         % this is on the NBS scale
        ./fH(F);                     % convert to SWS scale
end
F=(WhichKs==8);
if any(F)	
	% PURE WATER CASE
    % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
    % K1 from refit data from Harned and Davis,
    % J American Chemical Society, 65:2030-2037, 1943.
    % K2 from refit data from Harned and Scholes,
    % J American Chemical Society, 43:1706-1709, 1941.
	% This is only to be used for Sal=0 water (note the absence of S in the below formulations)
    % These are the thermodynamic Constants:
    lnK1(F) = 290.9097 - 14554.21./TempK(F) - 45.0575.*logTempK(F);
    K1(F) = exp(lnK1(F));
    lnK2(F) = 207.6548 - 11843.79./TempK(F) - 33.6485.*logTempK(F);
    K2(F) = exp(lnK2(F));
end
F=(WhichKs==9);
if any(F)
    % From Cai and Wang 1998, for estuarine use.
	% Data used in this work is from:
	% K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	% K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	% Sigma of residuals between fits and above data: ±0.015, +0.040 for K1 and K2, respectively.
	% Sal 0-40, Temp 0.2-30
  % Limnol. Oceanogr. 43(4) (1998) 657-668
	% On the NBS scale
	% Their check values for F1 don't work out, not sure if this was correctly published...
	F1 = 200.1./TempK(F) + 0.3220;
	pK1(F) = 3404.71./TempK(F) + 0.032786.*TempK(F) - 14.8435 - 0.071692.*F1.*Sal(F).^0.5 + 0.0021487.*Sal(F);
    K1(F)  = 10.^-pK1(F)...         % this is on the NBS scale
        ./fH(F);                    % convert to SWS scale (uncertain at low Sal due to junction potential);
	F2 = -129.24./TempK(F) + 1.4381;
	pK2(F) = 2902.39./TempK(F) + 0.02379.*TempK(F) - 6.4980 - 0.3191.*F2.*Sal(F).^0.5 + 0.0198.*Sal(F);
    K2(F)  = 10.^-pK2(F)...         % this is on the NBS scale
        ./fH(F);                    % convert to SWS scale (uncertain at low Sal due to junction potential); 
end
F=(WhichKs==10);
if any(F)
    % From Lueker, Dickson, Keeling, 2000
	% This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work. 
    % Mar. Chem. 70 (2000) 105-119
    % Total scale and kg-sw
    pK1(F) = 3633.86./TempK(F)-61.2172+9.6777.*log(TempK(F))-0.011555.*Sal(F)+0.0001152.*Sal(F).^2;
	K1(F)  = 10.^-pK1(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
    pK2(F) = 471.78./TempK(F)+25.929 -3.16967.*log(TempK(F))-0.01781 .*Sal(F)+0.0001122.*Sal(F).^2;
	K2(F)  = 10.^-pK2(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end
F=(WhichKs==11);
if any(F)
	% Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
	% sigma for pK1 is reported to be 0.0056
	% sigma for pK2 is reported to be 0.010
	% This is from the abstract and pages 2536-2537
    pK1 =  -43.6977 - 0.0129037.*Sal(F) + 1.364e-4.*Sal(F).^2 + 2885.378./TempK(F) +  7.045159.*log(TempK(F));
    pK2 = -452.0940 + 13.142162.*Sal(F) - 8.101e-4.*Sal(F).^2 + 21263.61./TempK(F) + 68.483143.*log(TempK(F))...
				+ (-581.4428.*Sal(F) + 0.259601.*Sal(F).^2)./TempK(F) - 1.967035.*Sal(F).*log(TempK(F));
	K1(F) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	K2(F) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==12);
if any(F)
	% Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
	% Calculated from overdetermined WOCE-era field measurements 
	% sigma for pK1 is reported to be 0.005
	% sigma for pK2 is reported to be 0.008
	% This is from page 1715
    pK1 =  6.359 - 0.00664.*Sal(F) - 0.01322.*TempC(F) + 4.989e-5.*TempC(F).^2;
    pK2 =  9.867 - 0.01314.*Sal(F) - 0.01904.*TempC(F) + 2.448e-5.*TempC(F).^2;
	K1(F) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	K2(F) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==13);
if any(F)
    % From Millero 2006 work on pK1 and pK2 from titrations
	% Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
    % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
	pK1_0 = -126.34048 + 6320.813./TempK(F) + 19.568224*log(TempK(F));
	A_1   = 13.4191*Sal(F).^0.5 + 0.0331.*Sal(F) - 5.33e-5.*Sal(F).^2;
	B_1   = -530.123*Sal(F).^0.5 - 6.103.*Sal(F);
	C_1   = -2.06950.*Sal(F).^0.5;
	pK1(F)= A_1 + B_1./TempK(F) + C_1.*log(TempK(F)) + pK1_0; % pK1 sigma = 0.0054
    K1(F) = 10.^-(pK1(F));
	pK2_0= -90.18333 + 5143.692./TempK(F) + 14.613358*log(TempK(F));	
	A_2   = 21.0894*Sal(F).^0.5 + 0.1248.*Sal(F) - 3.687e-4.*Sal(F).^2;
	B_2   = -772.483*Sal(F).^0.5 - 20.051.*Sal(F);
	C_2   = -3.3336.*Sal(F).^0.5;
	pK2(F)= A_2 + B_2./TempK(F) + C_2.*log(TempK(F)) + pK2_0; %pK2 sigma = 0.011
    K2(F) = 10.^-(pK2(F));
end
F=(WhichKs==14);
if any(F)
  % From Millero, 2010, also for estuarine use.
	% Marine and Freshwater Research, v. 61, p. 139-142.
	% Fits through compilation of real seawater titration results:
	% Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
	% Constants for K's on the SWS;
	% This is from page 141
	pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	% This is from their table 2, page 140.
	A1 = 13.4038.*Sal(F).^0.5 + 0.03206.*Sal(F) - 5.242e-5.*Sal(F).^2;
	B1 = -530.659.*Sal(F).^0.5 - 5.8210.*Sal(F);
	C1 = -2.0664*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1;
	% This is from page 141
	pK20 =  -90.18333 + 5143.692./TempK(F) + 14.613358.*log(TempK(F));
	% This is from their table 3, page 140.
	A2 = 21.3728.*Sal(F).^0.5 + 0.1218.*Sal(F) - 3.688e-4.*Sal(F).^2;
	B2 = -788.289.*Sal(F).^0.5 - 19.189.*Sal(F);
	C2 = -3.374.*Sal(F).^0.5;
	pK2 = pK20 + A2 + B2./TempK(F) + C2.*log(TempK(F));
	K2(F) = 10.^-pK2;
end
F=(WhichKs==15);
% Added by J. C. Orr on 4 Dec 2016
if any(F)
    % From Waters, Millero, and Woosley, 2014
	% Mar. Chem., 165, 66-67, 2014
  % Corrigendum to "The free proton concentration scale for seawater pH".
	% Effectively, this is an update of Millero (2010) formulation (WhichKs==14)
	% Constants for K's on the SWS;
	pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	A1 = 13.409160.*Sal(F).^0.5 + 0.031646.*Sal(F) - 5.1895e-5.*Sal(F).^2;
	B1 = -531.3642.*Sal(F).^0.5 - 5.713.*Sal(F);
	C1 = -2.0669166.*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1;
	pK20 =  -90.18333 + 5143.692./TempK(F) + 14.613358.*log(TempK(F));
	A2 = 21.225890.*Sal(F).^0.5 + 0.12450870.*Sal(F) - 3.7243e-4.*Sal(F).^2;
	B2 = -779.3444.*Sal(F).^0.5 - 19.91739.*Sal(F);
	C2 = -3.3534679.*Sal(F).^0.5;
	pK2 = pK20 + A2 + B2./TempK(F) + C2.*log(TempK(F));
	K2(F) = 10.^-pK2;
end
F=(WhichKs==16);
% Added by J. D. Sharp on 9 Jul 2020
if any(F)
    % From Sulpis et al, 2020
	% Ocean Science Discussions, 16, 847-862
    % This study uses overdeterminations of the carbonate system to
    % iteratively fit K1 and K2
    pK1(F) = 8510.63./TempK(F)-172.4493+26.32996.*log(TempK(F))-0.011555.*Sal(F)+0.0001152.*Sal(F).^2;
	K1(F)  = 10.^-pK1(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
    pK2(F) = 4226.23./TempK(F)-59.4636+9.60817.*log(TempK(F))-0.01781 .*Sal(F)+0.0001122.*Sal(F).^2;
	K2(F)  = 10.^-pK2(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end

F=(WhichKs==17);
% Added by J. D. Sharp on 15 Feb 2021
if any(F)
    % From Schockman & Byrne, 2021
	% Geochimica et Cosmochimica Acta, in press
    % This study uses spectrophotometric pH measurements to determine
    % K1*K2 with unprecedented precision, and presents a new
    % parameterization for K2 based on these determinations
    % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
    pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	A1 = 13.568513.*Sal(F).^0.5 + 0.031645.*Sal(F) - 5.3834e-5.*Sal(F).^2;
	B1 = -539.2304.*Sal(F).^0.5 - 5.635.*Sal(F);
	C1 = -2.0901396.*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1...               % this is on the total pH scale in mol/kg-sw
        ./SWStoTOT(F);                % convert to SWS pH scale
    % K2 is based on measurements of K1*K2:
    pK2 = 116.8067 - 3655.02./TempK(F) - 16.45817.*log(TempK(F)) + ...
        0.04523.*Sal(F) - 0.615.*Sal(F).^0.5 - 0.0002799.*Sal(F).^2 + ...
        4.969.*(Sal(F)./TempK(F));
    K2(F)  = 10.^-pK2...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end

%***************************************************************************
%CorrectKsForPressureNow:
% Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
%       the free scale) are on the SWS scale.
%       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
%       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
%       (the pH scales are the same in this case); the other Ks don't matter.
%
%
% No salinity dependence is given for the pressure coefficients here.
% It is assumed that the salinity is at or very near Sali = 35.
% These are valid for the SWS pH scale, but the difference between this and
% the total only yields a difference of .004 pH units at 1000 bars, much
% less than the uncertainties in the values.
%****************************************************************************
% The sources used are:
% Millero, 1995:
%       Millero, F. J., Thermodynamics of the carbon dioxide system in the
%       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
%       See table 9 and eqs. 90-92, p. 675.
%       TYPO: a factor of 10^3 was left out of the definition of Kappa
%       TYPO: the value of R given is incorrect with the wrong units
%       TYPO: the values of the a's for H2S and H2O are from the 1983
%                values for fresh water
%       TYPO: the value of a1 for B(OH)3 should be +.1622
%        Table 9 on p. 675 has no values for Si.
%       There are a variety of other typos in Table 9 on p. 675.
%       There are other typos in the paper, and most of the check values
%       given don't check.
% Millero, 1992:
%       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
%       CRC Press, 1992. See chapter 6.
%       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
%               79, and 96 have typos).
% Millero, 1983:
%       Millero, Frank J., Influence of pressure on chemical processes in
%       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
%       Chester, R., Academic Press, 1983.
%       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
%       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
%       these two are necessary to match the values given in Table 43.24
% Millero, 1979:
%       Millero, F. J., The thermodynamics of the carbon dioxide system
%       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
%       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
% Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
%       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
%       This matches the GEOSECS results and is in Edmond and Gieskes.
% Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
%       boric acid, and the pH of seawater, Limnology and Oceanography
%       13:403-417, 1968.
% Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
%       seawater with respect to calcium carbonate under in situ conditions,
%       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
%****************************************************************************
% These references often disagree and give different fits for the same thing.
% They are not always just an update either; that is, Millero, 1995 may agree
%       with Millero, 1979, but differ from Millero, 1983.
% For WhichKs% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
%       KP3, and KSi as for the other cases. Peng et al didn't consider the
%       case of P different from 0. GEOSECS did consider pressure, but didn't
%       include Phos, Si, or OH, so including the factors here won't matter.
% For WhichKs% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
%       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
%       including the factors won't matter.
%****************************************************************************
%       deltaVs are in cm3/mole
%       Kappas are in cm3/mole/bar
%****************************************************************************

%CorrectK1K2KBForPressure:
deltaV    = nan(ntps,1); Kappa     = nan(ntps,1);
lnK1fac   = nan(ntps,1); lnK2fac   = nan(ntps,1);
lnKBfac   = nan(ntps,1);
F=(WhichKs==8);
if any(F)
    %***PressureEffectsOnK1inFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  = -30.54 + 0.1849 .*TempC(F) - 0.0023366.*TempC(F).^2;
    Kappa(F)   = (-6.22 + 0.1368 .*TempC(F) - 0.001233 .*TempC(F).^2)./1000;
    lnK1fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %***PressureEffectsOnK2inFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  = -29.81 + 0.115.*TempC(F) - 0.001816.*TempC(F).^2;
    Kappa(F)   = (-5.74 + 0.093.*TempC(F) - 0.001896.*TempC(F).^2)./1000;
    lnK2fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    lnKBfac(F) = 0 ;%; this doesn't matter since TB = 0 for this case
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
    %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
    %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
    %               but the fits are the same as those in
    %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
    %               who in turn quote Li, personal communication
    lnK1fac(F) = (24.2 - 0.085.*TempC(F)).*Pbar(F)./RT(F);
    lnK2fac(F) = (16.4 - 0.04 .*TempC(F)).*Pbar(F)./RT(F);
    %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
    %               and matches the GEOSECS results
    lnKBfac(F) = (27.5 - 0.095.*TempC(F)).*Pbar(F)./RT(F);
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    %***PressureEffectsOnK1:
    %               These are from Millero, 1995.
    %               They are the same as Millero, 1979 and Millero, 1992.
    %               They are from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -25.5 + 0.1271.*TempC(F);
    %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = (-3.08 + 0.0877.*TempC(F))./1000;
    %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
 	lnK1fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               The fits given in Millero, 1983 are somewhat different.
    
    %***PressureEffectsOnK2:
    %               These are from Millero, 1995.
    %               They are the same as Millero, 1979 and Millero, 1992.
    %               They are from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -15.82 - 0.0219.*TempC(F);
    %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = (1.13 - 0.1475.*TempC(F))./1000;
    %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
	lnK2fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               The fit given in Millero, 1983 is different.
    %               Not by a lot for deltaV, but by much for Kappa. %
    
    %***PressureEffectsOnKB:
    %               This is from Millero, 1979.
    %               It is from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -29.48 + 0.1622.*TempC(F) - 0.002608.*TempC(F).^2;
    %               Millero, 1983 has:
    %                 'deltaV = -28.56 + .1211.*TempCi - .000321.*TempCi.*TempCi
    %               Millero, 1992 has:
    %                 'deltaV = -29.48 + .1622.*TempCi + .295.*(Sali - 34.8)
    %               Millero, 1995 has:
    %                 'deltaV = -29.48 - .1622.*TempCi - .002608.*TempCi.*TempCi
    %                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = -2.84./1000; % Millero, 1979
    %               Millero, 1992 and Millero, 1995 also have this.
    %                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
    %               Millero, 1983 has:
    %                 'Kappa = (-3 + .0427.*TempCi)./1000
    lnKBfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
end

% CorrectKWForPressure:
lnKWfac   = nan(ntps,1);
F=(WhichKs==8);
if any(F)
    % PressureEffectsOnKWinFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  =  -25.6 + 0.2324.*TempC(F) - 0.0036246.*TempC(F).^2;
    Kappa(F)   = (-7.33 + 0.1368.*TempC(F) - 0.001233 .*TempC(F).^2)./1000;
 	lnKWfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);

    %               NOTE the temperature dependence of KappaK1 and KappaKW
    %               for fresh water in Millero, 1983 are the same.
end
F=(WhichKs~=8);
if any(F)
    % GEOSECS doesn't include OH term, so this won't matter.
    % Peng et al didn't include pressure, but here I assume that the KW correction
    %       is the same as for the other seawater cases.
    % PressureEffectsOnKW:
    %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
    deltaV(F)  = -20.02 + 0.1119.*TempC(F) - 0.001409.*TempC(F).^2;
    %               Millero, 1992 and Millero, 1995 have:
    Kappa(F)   = (-5.13 + 0.0794.*TempC(F))./1000; % Millero, 1983
    %               Millero, 1995 has this too, but Millero, 1992 is different.
	lnKWfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               Millero, 1979 does not list values for these.
end

% PressureEffectsOnKF:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KF is on the free pH scale.
deltaV = -9.78 - 0.009.*TempC - 0.000942.*TempC.^2;
Kappa = (-3.91 + 0.054.*TempC)./1000;
lnKFfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKS:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KS is on the free pH scale.
deltaV = -18.03 + 0.0466.*TempC + 0.000316.*TempC.^2;
Kappa = (-4.53 + 0.09.*TempC)./1000;
lnKSfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKP1KP2KP3KSiForPressure:
% These corrections don't matter for the GEOSECS choice (WhichKs% = 6) and
%       the freshwater choice (WhichKs% = 8). For the Peng choice I assume
%       that they are the same as for the other choices (WhichKs% = 1 to 5).
% The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
%       same as Millero, 1983.
% PressureEffectsOnKP1:
deltaV = -14.51 + 0.1211.*TempC - 0.000321.*TempC.^2;
Kappa  = (-2.67 + 0.0427.*TempC)./1000;
lnKP1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP2:
deltaV = -23.12 + 0.1758.*TempC - 0.002647.*TempC.^2;
Kappa  = (-5.15 + 0.09  .*TempC)./1000;
lnKP2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP3:
deltaV = -26.57 + 0.202 .*TempC - 0.003042.*TempC.^2;
Kappa  = (-4.08 + 0.0714.*TempC)./1000;
lnKP3fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKSi:
%  The only mention of this is Millero, 1995 where it is stated that the
%    values have been estimated from the values of boric acid. HOWEVER,
%    there is no listing of the values in the table.
%    I used the values for boric acid from above.
deltaV = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;
Kappa  = -2.84./1000;
lnKSifac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKNH4: added by J. Sharp
deltaV = -26.43 + 0.0889.*TempC - 0.000905.*TempC.^2;
Kappa  = (-5.03 + 0.0814.*TempC)./1000;
lnKNH4fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKH2S: added by J. Sharp
deltaV = -11.07 - 0.009.*TempC - 0.000942.*TempC.^2;
Kappa  = (-2.89 + 0.054 .*TempC)./1000;
lnKH2Sfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKsForPressureHere:
K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
KWfac  = exp(lnKWfac);  KW  = KW .*KWfac;
KBfac  = exp(lnKBfac);  KB  = KB .*KBfac;
KFfac  = exp(lnKFfac);  KF  = KF .*KFfac;
KSfac  = exp(lnKSfac);  KS  = KS .*KSfac;
KP1fac = exp(lnKP1fac); KP1 = KP1.*KP1fac;
KP2fac = exp(lnKP2fac); KP2 = KP2.*KP2fac;
KP3fac = exp(lnKP3fac); KP3 = KP3.*KP3fac;
KSifac = exp(lnKSifac); KSi = KSi.*KSifac;
KNH4fac= exp(lnKNH4fac);KNH4= KNH4.*KNH4fac; % added by J. Sharp
KH2Sfac= exp(lnKH2Sfac);KH2S= KH2S.*KH2Sfac; % added by J. Sharp

% CorrectpHScaleConversionsForPressure:
% fH has been assumed to be independent of pressure.
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

%  The values KS and KF are already pressure-corrected, so the pH scale
%  conversions are now valid at pressure.

% FindpHScaleConversionFactor:
% this is the scale they will be put on
pHfactor  = nan(ntps,1);
F=(pHScale==1); %Total
pHfactor(F) = SWStoTOT(F);
F=(pHScale==2); %SWS, they are all on this now
pHfactor(F) = 1;
F=(pHScale==3); %pHfree
pHfactor(F) = SWStoTOT(F)./FREEtoTOT(F);
F=(pHScale==4); %pHNBS
pHfactor(F) = fH(F);

% ConvertFromSWSpHScaleToChosenScale:
K1   = K1.* pHfactor;  K2   = K2.* pHfactor;
KW   = KW.* pHfactor;  KB   = KB.* pHfactor;
KP1  = KP1.*pHfactor;  KP2  = KP2.*pHfactor;
KP3  = KP3.*pHfactor;  KSi  = KSi.*pHfactor;
KNH4 = KNH4.*pHfactor; KH2S = KH2S.*pHfactor;

% CalculateFugacityConstants:
% In previos versions of CO2SYS, the fugacity factor was calculated
% assuming pressure at one atmosphere, or close to it. Starting with
% v3.2.1, an option to use in situ pressure is provided.
%       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
%       Delta and B in cm3/mol
FugFac=ones(ntps,1);
Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
% For a mixture of CO2 and air at in situ pressure;
xc2 = 1; % assumed to be 1, though not strictly correct (xc2 = [1-xCO2]^2)
P1atm = 1.01325; % atmospheric pressure in bar
if p_opt == 0
    FugFac = exp((b + 2.*xc2.*Delta).*P1atm./RT); % FugFac at 1 atm
elseif p_opt == 1
    FugFac = exp((b + 2.*xc2.*Delta).*(P1atm+Pbar)./RT); % FugFac at in situ pressure
else
    disp('co2_press must be set to either 0 or 1'); % Display error message
end
F=(WhichKs==6 | WhichKs==7); % GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
FugFac(F) = 1;

% CalculateVPFac:
% Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
%       seawater, Marine Chemistry 8:347-359, 1980.
% They fit the data of Goff and Gratch (1946) with the vapor pressure
%       lowering by sea salt as given by Robinson (1954).
% This fits the more complicated Goff and Gratch, and Robinson equations
%       from 273 to 313 deg K and 0 to 40 Sali with a standard error
%       of .015%, about 5 uatm over this range.
% This may be on IPTS-29 since they didn't mention the temperature scale,
%       and the data of Goff and Gratch came before IPTS-48.
% The references are:
% Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
%       to 212 deg F, Transactions of the American Society of Heating and
%       Ventilating Engineers 52:95-122, 1946.
% Robinson, Journal of the Marine Biological Association of the U. K.
%       33:449-455, 1954.
%       This is eq. 10 on p. 350.
%       This is in atmospheres.
VPWP = exp(24.4543 - 67.4509.*(100./TempK) - 4.8489.*log(TempK./100));
VPCorrWP = exp(-0.000544.*Sal);
VPSWWP = VPWP.*VPCorrWP;
VPFac = 1 - VPSWWP; % this assumes 1 atmosphere
end % end nested function