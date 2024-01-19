
function Ks = calculate_equilibrium_constants(TempC,Pdbar,pH_scale,p_opt,gas_constant)
    global which_k1_k2_constants_GLOBAL which_kso4_constant_GLOBAL which_kf_constant_GLOBAL which_boron_GLOBAL Pbar;
    global fH ntps temp_k_GLOBAL log_temp_k_GLOBAL;
    global boron_concentration_GLOBAL fluorine_concentration_GLOBAL sulphate_concentration_GLOBAL CAL salinity_GLOBAL;
    
    % SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
    % Inputs: pHScale%, which_k1_k2_constants_GLOBAL%, which_kso4_constant_GLOBAL%, Sali, temperature_in_GLOBAL, Pdbar
    % Outputs: K0, K(), T(), fH
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
    temp_k_GLOBAL    = TempC + 273.15;
    log_temp_k_GLOBAL = log(temp_k_GLOBAL);
    Pbar     = Pdbar ./ 10;
    
    % Generate empty vectors for holding results
    boron_concentration_GLOBAL = nan(ntps,1);
    fluorine_concentration_GLOBAL = nan(ntps,1);
    sulphate_concentration_GLOBAL = nan(ntps,1);
    CAL = nan(ntps,1);
    
    % CalculateTB - Total Borate:
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8); % Pure water case.
    if any(selected_GLOBAL)
        boron_concentration_GLOBAL(selected_GLOBAL) = 0;
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        boron_concentration_GLOBAL(selected_GLOBAL) = 0.0004106.*salinity_GLOBAL(selected_GLOBAL)./35; % in mol/kg-SW
        % this is .00001173.*Sali
        % this is about 1% lower than Uppstrom's value
        % Culkin, selected_GLOBAL., in Chemical Oceanography,
        % ed. Riley and Skirrow, 1965:
        % GEOSECS references this, but this value is not explicitly
        % given here
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8); % All other cases
    if any(selected_GLOBAL)
	    FF=selected_GLOBAL&(which_boron_GLOBAL==1); % If user opted for Uppstrom's values:
	    if any(FF)
	        % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
	        % this is .000416.*Sali./35. = .0000119.*Sali
		    % boron_concentration_GLOBAL(FF) = (0.000232./10.811).*(salinity_GLOBAL(FF)./1.80655); % in mol/kg-SW
	        boron_concentration_GLOBAL(FF) =  0.0004157.*salinity_GLOBAL(FF)./35; % in mol/kg-SW
	    end
	    FF=selected_GLOBAL&(which_boron_GLOBAL==2); % If user opted for the Lee et al. values:
	    if any(FF)
		    % Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.	
	 	    % Geochimica Et Cosmochimica Acta 74 (6): 1801-1811.
		    boron_concentration_GLOBAL(FF) =  0.0004326.*salinity_GLOBAL(FF)./35; % in mol/kg-SW
	    end
    end
    
    % CalculateCAL - Total Calcium:
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7);
        % Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
        % this is .010285.*Sali./35
        CAL(selected_GLOBAL) = 0.02128./40.087.*(salinity_GLOBAL(selected_GLOBAL)./1.80655);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
        % *** CalculateCaforGEOSECS:
        % Culkin, selected_GLOBAL, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982) 
        CAL(selected_GLOBAL) = 0.01026.*salinity_GLOBAL(selected_GLOBAL)./35;
    
    % CalculateTF;
    % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
    % this is .000068.*Sali./35. = .00000195.*Sali
    fluorine_concentration_GLOBAL = (0.000067./18.998).*(salinity_GLOBAL./1.80655); % in mol/kg-SW
    
    % CalculateTS ;
    % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
    % this is .02824.*Sali./35. = .0008067.*Sali
    sulphate_concentration_GLOBAL = (0.14./96.062).*(salinity_GLOBAL./1.80655); % in mol/kg-SW
    
    % CalculateK0:
    % Weiss, R. selected_GLOBAL., Marine Chemistry 2:203-215, 1974.
    TempK100  = temp_k_GLOBAL./100;
    lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + salinity_GLOBAL .*...
        (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
    K0   = exp(lnK0);                 % this is in mol/kg-SW/atm
    vCO2 = 32.3;                      % partial molal volume of CO2 (cm3 / mol)
                                      % from Weiss (1974, Appendix, paragraph 3)
    if p_opt == 0
        prscorr = 1; % Set pressure correction to 1
    elseif p_opt == 1
        prscorr = exp((-Pbar).*vCO2./(gas_constant.*temp_k_GLOBAL)); % Calculate pressure correction to K0
    else
        disp('co2_press must be set to either 0 or 1'); % Display error message
    end         
    K0   = K0 .* prscorr; % this is in mol/kg-SW/atm
    
    % CalculateIonS:
    % This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
    IonS         = 19.924 .* salinity_GLOBAL ./ (1000 - 1.005   .* salinity_GLOBAL);
    
    % CalculateKS:
    lnKS   = nan(ntps,1); pKS  = nan(ntps,1); KS   = nan(ntps,1);
    logKS0 = nan(ntps,1); logKSK0 = nan(ntps,1);
    selected_GLOBAL=(which_kso4_constant_GLOBAL==1);
    if any(selected_GLOBAL)
        % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
        % The goodness of fit is .021.
        % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
        % TYPO on p. 121: the constant e9 should be e8.
        % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
      lnKS(selected_GLOBAL) = -4276.1./temp_k_GLOBAL(selected_GLOBAL) + 141.328 - 23.093.*log_temp_k_GLOBAL(selected_GLOBAL) +...             
          (-13856./temp_k_GLOBAL(selected_GLOBAL) + 324.57 - 47.986.*log_temp_k_GLOBAL(selected_GLOBAL)).*sqrt(IonS(selected_GLOBAL)) +...     
          (35474./temp_k_GLOBAL(selected_GLOBAL) - 771.54 + 114.723.*log_temp_k_GLOBAL(selected_GLOBAL)).*IonS(selected_GLOBAL) +...           
          (-2698./temp_k_GLOBAL(selected_GLOBAL)).*sqrt(IonS(selected_GLOBAL)).*IonS(selected_GLOBAL) + (1776./temp_k_GLOBAL(selected_GLOBAL)).*IonS(selected_GLOBAL).^2; 
	    KS(selected_GLOBAL) = exp(lnKS(selected_GLOBAL))...            % this is on the free pH scale in mol/kg-H2O
            .* (1 - 0.001005 .* salinity_GLOBAL(selected_GLOBAL));   % convert to mol/kg-SW
    end
    selected_GLOBAL=(which_kso4_constant_GLOBAL==2);
    if any(selected_GLOBAL)
        % Khoo et al, Analytical Chemistry, 49(1):29-34, 1977
        % KS was found by titrations with a hydrogen electrode
        % of artificial seawater containing sulfate (but without selected_GLOBAL)
        % at 3 salinities from 20 to 45 and artificial seawater NOT
        % containing sulfate (nor selected_GLOBAL) at 16 salinities from 15 to 45,
        % both at temperatures from 5 to 40 deg C.
        % KS is on the Free pH scale (inherently so).
        % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
        % He finds log(beta) which = my pKS;
        % his beta is an association constant.
        % The rms error is .0021 in pKS, or about .5% in KS.
        % This is equation 20 on p. 33:
        pKS(selected_GLOBAL) = 647.59 ./ temp_k_GLOBAL(selected_GLOBAL) - 6.3451 + 0.019085.*temp_k_GLOBAL(selected_GLOBAL) - 0.5208.*sqrt(IonS(selected_GLOBAL));
        KS(selected_GLOBAL) = 10.^(-pKS(selected_GLOBAL))...          % this is on the free pH scale in mol/kg-H2O
            .* (1 - 0.001005.*salinity_GLOBAL(selected_GLOBAL));    % convert to mol/kg-SW
    end
    selected_GLOBAL=(which_kso4_constant_GLOBAL==3);
    if any(selected_GLOBAL)
        % Waters and Millero, Marine Chemistry, 149: 8-22, 2013, with corrections from
        % Waters et al, Marine Chemistry, 165: 66-67, 2014
        logKS0(selected_GLOBAL) = 562.69486 - 102.5154.*log_temp_k_GLOBAL(selected_GLOBAL) - 0.0001117033.*temp_k_GLOBAL(selected_GLOBAL).*temp_k_GLOBAL(selected_GLOBAL) + ...
            0.2477538.*temp_k_GLOBAL(selected_GLOBAL) - 13273.76./temp_k_GLOBAL(selected_GLOBAL);
        logKSK0(selected_GLOBAL) = (4.24666 - 0.152671.*temp_k_GLOBAL(selected_GLOBAL) + 0.0267059.*temp_k_GLOBAL(selected_GLOBAL).*log_temp_k_GLOBAL(selected_GLOBAL) - 0.000042128.*temp_k_GLOBAL(selected_GLOBAL).*temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^0.5 + ...
            (0.2542181 - 0.00509534.*temp_k_GLOBAL(selected_GLOBAL) + 0.00071589.*temp_k_GLOBAL(selected_GLOBAL).*log_temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL) + (-0.00291179 + 0.0000209968.*temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^1.5 + ...
            -0.0000403724.*salinity_GLOBAL(selected_GLOBAL).^2;
        KS(selected_GLOBAL) = ((10.^(logKSK0(selected_GLOBAL))).*(10.^logKS0(selected_GLOBAL))) ... % this is on the free pH scale in mol/kg-H2O
            .* (1 - 0.001005.*salinity_GLOBAL(selected_GLOBAL));                    % convert to mol/kg-SW
    end
    
    % CalculateKF:
    KF = NaN(ntps, 1);  % added preallocation here and selected_GLOBAL-indexing below // MPH
    selected_GLOBAL=(which_kf_constant_GLOBAL==1);
    if any(selected_GLOBAL)
        % Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
        lnKF = 1590.2./temp_k_GLOBAL - 12.641 + 1.525.*IonS.^0.5;
        KF(selected_GLOBAL)   = exp(lnKF(selected_GLOBAL))...                 % this is on the free pH scale in mol/kg-H2O
            .*(1 - 0.001005.*salinity_GLOBAL(selected_GLOBAL));          % convert to mol/kg-SW
    end
    selected_GLOBAL=(which_kf_constant_GLOBAL==2);
    if any(selected_GLOBAL)
        % Perez and Fraga 1987 (to be used for S: 10-40, T: 9-33)
        % P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
        lnKF = 874./temp_k_GLOBAL - 9.68 + 0.111.*salinity_GLOBAL.^0.5;
        KF(selected_GLOBAL)   = exp(lnKF(selected_GLOBAL));                   % this is on the free pH scale in mol/kg-SW
    end
    
    % CalculatepHScaleConversionFactors:
    %       These are NOT pressure-corrected.
    SWStoTOT  = (1 + sulphate_concentration_GLOBAL./KS)./(1 + sulphate_concentration_GLOBAL./KS + fluorine_concentration_GLOBAL./KF);
    FREEtoTOT =  1 + sulphate_concentration_GLOBAL./KS;
    
    % CalculatefH
    fH = nan(ntps,1);
    % Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8);
    if any(selected_GLOBAL)
        fH(selected_GLOBAL) = 1; % this shouldn't occur in the program for this case
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        fH(selected_GLOBAL) = 1.29 - 0.00204.*  temp_k_GLOBAL(selected_GLOBAL) + (0.00046 -...
            0.00000148.*temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).*salinity_GLOBAL(selected_GLOBAL);
        % Peng et al, Tellus 39B:439-458, 1987:
        % They reference the GEOSECS report, but round the value
        % given there off so that it is about .008 (1%) lower. It
        % doesn't agree with the check value they give on p. 456.
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8);
    if any(selected_GLOBAL)
        fH(selected_GLOBAL) = 1.2948 - 0.002036.*temp_k_GLOBAL(selected_GLOBAL) + (0.0004607 -...
            0.000001475.*temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^2;
        % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
        % v. 3, 1982 (p. 80);
    end
    
    % CalculateKB:
    KB      = nan(ntps,1); logKB   = nan(ntps,1);
    lnKBtop = nan(ntps,1); lnKB    = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8); % Pure water case
    if any(selected_GLOBAL)
        KB(selected_GLOBAL) = 0;
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        % This is for GEOSECS and Peng et al.
        % Lyman, John, UCLA Thesis, 1957
        % fit by Li et al, JGR 74:5507-5525, 1969:
        logKB(selected_GLOBAL) = -9.26 + 0.00886.*salinity_GLOBAL(selected_GLOBAL) + 0.01.*TempC(selected_GLOBAL);
        KB(selected_GLOBAL) = 10.^(logKB(selected_GLOBAL))...  % this is on the NBS scale
            ./fH(selected_GLOBAL);               % convert to the SWS scale
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8);
    if any(selected_GLOBAL)
        % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
        lnKBtop(selected_GLOBAL) = -8966.9 - 2890.53.*sqrt(salinity_GLOBAL(selected_GLOBAL)) - 77.942.*salinity_GLOBAL(selected_GLOBAL) +...
            1.728.*sqrt(salinity_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL) - 0.0996.*salinity_GLOBAL(selected_GLOBAL).^2;
        lnKB(selected_GLOBAL) = lnKBtop(selected_GLOBAL)./temp_k_GLOBAL(selected_GLOBAL) + 148.0248 + 137.1942.*sqrt(salinity_GLOBAL(selected_GLOBAL)) +...
            1.62142.*salinity_GLOBAL(selected_GLOBAL) + (-24.4344 - 25.085.*sqrt(salinity_GLOBAL(selected_GLOBAL)) - 0.2474.*...
            salinity_GLOBAL(selected_GLOBAL)).*log_temp_k_GLOBAL(selected_GLOBAL) + 0.053105.*sqrt(salinity_GLOBAL(selected_GLOBAL)).*temp_k_GLOBAL(selected_GLOBAL);
        KB(selected_GLOBAL) = exp(lnKB(selected_GLOBAL))...    % this is on the total pH scale in mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);         % convert to SWS pH scale
    end
    
    % CalculateKW:
    lnKW = nan(ntps,1); KW = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
        lnKW(selected_GLOBAL) = 148.9802 - 13847.26./temp_k_GLOBAL(selected_GLOBAL) - 23.6521.*log_temp_k_GLOBAL(selected_GLOBAL) +...
            (-79.2447 + 3298.72./temp_k_GLOBAL(selected_GLOBAL) + 12.0408.*log_temp_k_GLOBAL(selected_GLOBAL)).*...
            sqrt(salinity_GLOBAL(selected_GLOBAL)) - 0.019813.*salinity_GLOBAL(selected_GLOBAL);
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8);
    if any(selected_GLOBAL)
        % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
        % refit data of Harned and Owen, The Physical Chemistry of
        % Electrolyte Solutions, 1958
        lnKW(selected_GLOBAL) = 148.9802 - 13847.26./temp_k_GLOBAL(selected_GLOBAL) - 23.6521.*log_temp_k_GLOBAL(selected_GLOBAL);
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8);
    if any(selected_GLOBAL)
        % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
        % his check value of 1.6 umol/kg-SW should be 6.2
        lnKW(selected_GLOBAL) = 148.9802 - 13847.26./temp_k_GLOBAL(selected_GLOBAL) - 23.6521.*log_temp_k_GLOBAL(selected_GLOBAL) +...
            (-5.977 + 118.67./temp_k_GLOBAL(selected_GLOBAL) + 1.0495.*log_temp_k_GLOBAL(selected_GLOBAL)).*...
            sqrt(salinity_GLOBAL(selected_GLOBAL)) - 0.01615.*salinity_GLOBAL(selected_GLOBAL);
    end
    KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6);
    if any(selected_GLOBAL)
        KW(selected_GLOBAL) = 0; % GEOSECS doesn't include OH effects
    end
    
    % CalculateKP1KP2KP3KSi:
    KP1      = nan(ntps,1); KP2      = nan(ntps,1);
    KP3      = nan(ntps,1); KSi      = nan(ntps,1);
    lnKP1    = nan(ntps,1); lnKP2    = nan(ntps,1);
    lnKP3    = nan(ntps,1); lnKSi    = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        KP1(selected_GLOBAL) = 0.02;
        % Peng et al don't include the contribution from this term,
        % but it is so small it doesn't contribute. It needs to be
        % kept so that the routines work ok.
        % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
        % Limnology and Oceanography 12:243-252, 1967:
        % these are only for sals 33 to 36 and are on the NBS scale
        KP2(selected_GLOBAL) = exp(-9.039 - 1450./temp_k_GLOBAL(selected_GLOBAL))... % this is on the NBS scale
            ./fH(selected_GLOBAL);                          % convert to SWS scale
        KP3(selected_GLOBAL) = exp(4.466 - 7276./temp_k_GLOBAL(selected_GLOBAL))...  % this is on the NBS scale
            ./fH(selected_GLOBAL);                          % convert to SWS scale
        % Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
        % The Chemical Society (London), Special Publ. 17:751, 1964:
        KSi(selected_GLOBAL) = 0.0000000004...              % this is on the NBS scale
            ./fH(selected_GLOBAL);                          % convert to SWS scale
        
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==8);
    if any(selected_GLOBAL)
        KP1(selected_GLOBAL) = 0; KP2(selected_GLOBAL) = 0; KP3(selected_GLOBAL) = 0; KSi(selected_GLOBAL) = 0;
        % Neither the GEOSECS choice nor the freshwater choice
        % include contributions from phosphate or silicate.
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8);
    if any(selected_GLOBAL)
        % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
        % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
        % KSi was given on the SWS pH scale in molal units.
        lnKP1(selected_GLOBAL) = -4576.752./temp_k_GLOBAL(selected_GLOBAL) + 115.54 - 18.453.*log_temp_k_GLOBAL(selected_GLOBAL) + (-106.736./temp_k_GLOBAL(selected_GLOBAL) +...
            0.69171).*sqrt(salinity_GLOBAL(selected_GLOBAL)) + (-0.65643./temp_k_GLOBAL(selected_GLOBAL) - 0.01844).*salinity_GLOBAL(selected_GLOBAL);
        KP1(selected_GLOBAL) = exp(lnKP1(selected_GLOBAL));
        lnKP2(selected_GLOBAL) = -8814.715./temp_k_GLOBAL(selected_GLOBAL) + 172.1033 - 27.927.*log_temp_k_GLOBAL(selected_GLOBAL) + (-160.34./temp_k_GLOBAL(selected_GLOBAL) +...
            1.3566).*sqrt(salinity_GLOBAL(selected_GLOBAL)) + (0.37335./temp_k_GLOBAL(selected_GLOBAL) - 0.05778).*salinity_GLOBAL(selected_GLOBAL);
        KP2(selected_GLOBAL) = exp(lnKP2(selected_GLOBAL));
        lnKP3(selected_GLOBAL) = -3070.75./temp_k_GLOBAL(selected_GLOBAL) - 18.126 + (17.27039./temp_k_GLOBAL(selected_GLOBAL) + 2.81197).*sqrt(salinity_GLOBAL(selected_GLOBAL)) +...
            (-44.99486./temp_k_GLOBAL(selected_GLOBAL) - 0.09984).*salinity_GLOBAL(selected_GLOBAL);
        KP3(selected_GLOBAL) = exp(lnKP3(selected_GLOBAL));
        lnKSi(selected_GLOBAL) = -8904.2./temp_k_GLOBAL(selected_GLOBAL) + 117.4 - 19.334.*log_temp_k_GLOBAL(selected_GLOBAL) + (-458.79./temp_k_GLOBAL(selected_GLOBAL) +...
            3.5913).*sqrt(IonS(selected_GLOBAL)) + (188.74./temp_k_GLOBAL(selected_GLOBAL) - 1.5998).*IonS(selected_GLOBAL) +...
            (-12.1652./temp_k_GLOBAL(selected_GLOBAL) + 0.07871).*IonS(selected_GLOBAL).^2;
        KSi(selected_GLOBAL) = exp(lnKSi(selected_GLOBAL))...                % this is on the SWS pH scale in mol/kg-H2O
            .*(1 - 0.001005.*salinity_GLOBAL(selected_GLOBAL));        % convert to mol/kg-SW
    end
    
    % Calculate KNH4 and KH2S: added by J. Sharp
    KNH4           = nan(ntps,1); KH2S       = nan(ntps,1);
    PKNH4expCW     = nan(ntps,1); lnKH2S     = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7 | which_k1_k2_constants_GLOBAL==8); % GEOSECS or freshwater cases
    if any(selected_GLOBAL)
        KNH4(selected_GLOBAL) = 0;
        KH2S(selected_GLOBAL) = 0;
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8); % All other cases
    if any(selected_GLOBAL)
    % Ammonia dissociation constant from Yao and Millero (1995)
    %   KNH4(selected_GLOBAL) = (exp(-6285.33./temp_k_GLOBAL(selected_GLOBAL)+0.0001635.*temp_k_GLOBAL(selected_GLOBAL)-0.25444+...
    %             (0.46532-123.7184./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^0.5+(-0.01992+...
    %             3.17556./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL)))...
    %             ./SWStoTOT(selected_GLOBAL);                    % convert to SWS pH scale
    % Ammonia dissociation constant from Clegg and Whitfield (1995)
      PKNH4(selected_GLOBAL) = 9.244605-2729.33.*(1./298.15-1./temp_k_GLOBAL(selected_GLOBAL)) +...
              (0.04203362-11.24742./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^0.25+...  % added missing (selected_GLOBAL) index on salinity_GLOBAL // MPH
              (-13.6416+1.176949.*temp_k_GLOBAL(selected_GLOBAL).^0.5-...
              0.02860785.*temp_k_GLOBAL(selected_GLOBAL)+545.4834./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^0.5+...
              (-0.1462507+0.0090226468.*temp_k_GLOBAL(selected_GLOBAL).^0.5-...
              0.0001471361.*temp_k_GLOBAL(selected_GLOBAL)+10.5425./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^1.5+...
              (0.004669309-0.0001691742.*temp_k_GLOBAL(selected_GLOBAL).^0.5-...
              0.5677934./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^2+...
              (-2.354039E-05+0.009698623./temp_k_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL).^2.5;
      KNH4(selected_GLOBAL)  = 10.^-PKNH4(selected_GLOBAL);                    % total scale, mol/kg-H2O
      KNH4(selected_GLOBAL)  = KNH4(selected_GLOBAL).*(1-0.001005.*salinity_GLOBAL(selected_GLOBAL)); % mol/kg-SW
      KNH4(selected_GLOBAL)  = KNH4(selected_GLOBAL)./SWStoTOT(selected_GLOBAL);             % converts to SWS pH scale
    
    % First hydrogen sulfide dissociation constant from Millero et al. (1988)
      KH2S(selected_GLOBAL)  = (exp(225.838-13275.3./temp_k_GLOBAL(selected_GLOBAL)-34.6435.*log(temp_k_GLOBAL(selected_GLOBAL))+...
                  0.3449.*salinity_GLOBAL(selected_GLOBAL).^0.5-0.0274.*salinity_GLOBAL(selected_GLOBAL)))...
                  ./SWStoTOT(selected_GLOBAL);                    % convert to SWS pH scale
    
    end
    
    % CalculateK1K2:
    logK1    = nan(ntps,1); lnK1     = nan(ntps,1);
    pK1      = nan(ntps,1); K1       = nan(ntps,1);
    logK2    = nan(ntps,1); lnK2     = nan(ntps,1);
    pK2      = nan(ntps,1); K2       = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==1);
    if any(selected_GLOBAL)
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
        lnK1(selected_GLOBAL) = 2.83655 - 2307.1266./temp_k_GLOBAL(selected_GLOBAL) - 1.5529413.*log_temp_k_GLOBAL(selected_GLOBAL) +...
            (-0.20760841 - 4.0484./temp_k_GLOBAL(selected_GLOBAL)).*sqrt(salinity_GLOBAL(selected_GLOBAL)) + 0.08468345.*salinity_GLOBAL(selected_GLOBAL) -...
            0.00654208.*sqrt(salinity_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL);
        K1(selected_GLOBAL) = exp(lnK1(selected_GLOBAL))...            % this is on the total pH scale in mol/kg-H2O
            .*(1 - 0.001005.*salinity_GLOBAL(selected_GLOBAL))...    % convert to mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                 % convert to SWS pH scale
        % This is eq. 30 on p. 254 and what they use in their abstract:
        lnK2(selected_GLOBAL) = -9.226508 - 3351.6106./temp_k_GLOBAL(selected_GLOBAL) - 0.2005743.*log_temp_k_GLOBAL(selected_GLOBAL) +...
            (-0.106901773 - 23.9722./temp_k_GLOBAL(selected_GLOBAL)).*sqrt(salinity_GLOBAL(selected_GLOBAL)) + 0.1130822.*salinity_GLOBAL(selected_GLOBAL) -...
            0.00846934.*sqrt(salinity_GLOBAL(selected_GLOBAL)).*salinity_GLOBAL(selected_GLOBAL);
        K2(selected_GLOBAL) = exp(lnK2(selected_GLOBAL))...            % this is on the total pH scale in mol/kg-H2O
            .*(1 - 0.001005.*salinity_GLOBAL(selected_GLOBAL))...    % convert to mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                 % convert to SWS pH scale
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==2);
    if any(selected_GLOBAL)
        % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
        % The 2s precision in pK1 is .011, or 2.5% in K1.
        % The 2s precision in pK2 is .02, or 4.5% in K2.
        % This is in Table 5 on p. 1652 and what they use in the abstract:
        pK1(selected_GLOBAL) = 812.27./temp_k_GLOBAL(selected_GLOBAL) + 3.356 - 0.00171.*salinity_GLOBAL(selected_GLOBAL).*log_temp_k_GLOBAL(selected_GLOBAL)...
            + 0.000091.*salinity_GLOBAL(selected_GLOBAL).^2;
        K1(selected_GLOBAL) = 10.^(-pK1(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
        % This is in Table 5 on p. 1652 and what they use in the abstract:
        pK2(selected_GLOBAL) = 1450.87./temp_k_GLOBAL(selected_GLOBAL) + 4.604 - 0.00385.*salinity_GLOBAL(selected_GLOBAL).*log_temp_k_GLOBAL(selected_GLOBAL)...
            + 0.000182.*salinity_GLOBAL(selected_GLOBAL).^2;
        K2(selected_GLOBAL) = 10.^(-pK2(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==3);
    if any(selected_GLOBAL)
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
        pK1(selected_GLOBAL) = 851.4./temp_k_GLOBAL(selected_GLOBAL) + 3.237 - 0.0106.*salinity_GLOBAL(selected_GLOBAL) + 0.000105.*salinity_GLOBAL(selected_GLOBAL).^2;
        K1(selected_GLOBAL) = 10.^(-pK1(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
        % This is from Table 4 on p. 1739.
        pK2(selected_GLOBAL) = -3885.4./temp_k_GLOBAL(selected_GLOBAL) + 125.844 - 18.141.*log_temp_k_GLOBAL(selected_GLOBAL)...
            - 0.0192.*salinity_GLOBAL(selected_GLOBAL) + 0.000132.*salinity_GLOBAL(selected_GLOBAL).^2;
        K2(selected_GLOBAL) = 10.^(-pK2(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==4);
    if any(selected_GLOBAL)
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
        pK1(selected_GLOBAL) = 3670.7./temp_k_GLOBAL(selected_GLOBAL) - 62.008 + 9.7944.*log_temp_k_GLOBAL(selected_GLOBAL)...
                 - 0.0118.*salinity_GLOBAL(selected_GLOBAL) + 0.000116.*salinity_GLOBAL(selected_GLOBAL).^2;
        K1(selected_GLOBAL) = 10.^(-pK1(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
        % This is in Table 4 on p. 1739.
        pK2(selected_GLOBAL) = 1394.7./temp_k_GLOBAL(selected_GLOBAL) + 4.777 - 0.0184.*salinity_GLOBAL(selected_GLOBAL) + 0.000118.*salinity_GLOBAL(selected_GLOBAL).^2;
        K2(selected_GLOBAL) = 10.^(-pK2(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==5);
    if any(selected_GLOBAL)
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
        pK1(selected_GLOBAL) = 845./temp_k_GLOBAL(selected_GLOBAL) + 3.248 - 0.0098.*salinity_GLOBAL(selected_GLOBAL) + 0.000087.*salinity_GLOBAL(selected_GLOBAL).^2;
        K1(selected_GLOBAL) = 10.^(-pK1(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
        % This is in Table 5 on p. 1740.
        pK2(selected_GLOBAL) = 1377.3./temp_k_GLOBAL(selected_GLOBAL) + 4.824 - 0.0185.*salinity_GLOBAL(selected_GLOBAL) + 0.000122.*salinity_GLOBAL(selected_GLOBAL).^2;
        K2(selected_GLOBAL) = 10.^(-pK2(selected_GLOBAL)); % this is on the SWS pH scale in mol/kg-SW
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
        % Limnology and Oceanography, 18(6):897-907, 1973.
	    % I.e., these are the original Mehrbach dissociation constants.
        % The 2s precision in pK1 is .005, or 1.2% in K1.
        % The 2s precision in pK2 is .008, or 2% in K2.
        pK1(selected_GLOBAL) = - 13.7201 + 0.031334.*temp_k_GLOBAL(selected_GLOBAL) + 3235.76./temp_k_GLOBAL(selected_GLOBAL)...
            + 1.3e-5*salinity_GLOBAL(selected_GLOBAL).*temp_k_GLOBAL(selected_GLOBAL) - 0.1032.*salinity_GLOBAL(selected_GLOBAL).^0.5;
        K1(selected_GLOBAL) = 10.^(-pK1(selected_GLOBAL))...         % this is on the NBS scale
            ./fH(selected_GLOBAL);                     % convert to SWS scale
        pK2(selected_GLOBAL) = 5371.9645 + 1.671221.*temp_k_GLOBAL(selected_GLOBAL) + 0.22913.*salinity_GLOBAL(selected_GLOBAL) + 18.3802.*log10(salinity_GLOBAL(selected_GLOBAL))...
                 - 128375.28./temp_k_GLOBAL(selected_GLOBAL) - 2194.3055.*log10(temp_k_GLOBAL(selected_GLOBAL)) - 8.0944e-4.*salinity_GLOBAL(selected_GLOBAL).*temp_k_GLOBAL(selected_GLOBAL)...
                 - 5617.11.*log10(salinity_GLOBAL(selected_GLOBAL))./temp_k_GLOBAL(selected_GLOBAL) + 2.136.*salinity_GLOBAL(selected_GLOBAL)./temp_k_GLOBAL(selected_GLOBAL); % pK2 is not defined for salinity_GLOBAL=0, since log10(0)=-inf
        K2(selected_GLOBAL) = 10.^(-pK2(selected_GLOBAL))...         % this is on the NBS scale
            ./fH(selected_GLOBAL);                     % convert to SWS scale
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8);
    if any(selected_GLOBAL)	
	    % PURE WATER CASE
        % Millero, selected_GLOBAL. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
        % K1 from refit data from Harned and Davis,
        % J American Chemical Society, 65:2030-2037, 1943.
        % K2 from refit data from Harned and Scholes,
        % J American Chemical Society, 43:1706-1709, 1941.
	    % This is only to be used for salinity_GLOBAL=0 water (note the absence of S in the below formulations)
        % These are the thermodynamic Constants:
        lnK1(selected_GLOBAL) = 290.9097 - 14554.21./temp_k_GLOBAL(selected_GLOBAL) - 45.0575.*log_temp_k_GLOBAL(selected_GLOBAL);
        K1(selected_GLOBAL) = exp(lnK1(selected_GLOBAL));
        lnK2(selected_GLOBAL) = 207.6548 - 11843.79./temp_k_GLOBAL(selected_GLOBAL) - 33.6485.*log_temp_k_GLOBAL(selected_GLOBAL);
        K2(selected_GLOBAL) = exp(lnK2(selected_GLOBAL));
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==9);
    if any(selected_GLOBAL)
        % From Cai and Wang 1998, for estuarine use.
	    % Data used in this work is from:
	    % K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	    % K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	    % Sigma of residuals between fits and above data: ±0.015, +0.040 for K1 and K2, respectively.
	    % salinity_GLOBAL 0-40, Temp 0.2-30
      % Limnol. Oceanogr. 43(4) (1998) 657-668
	    % On the NBS scale
	    % Their check values for F1 don't work out, not sure if this was correctly published...
	    F1 = 200.1./temp_k_GLOBAL(selected_GLOBAL) + 0.3220;
	    pK1(selected_GLOBAL) = 3404.71./temp_k_GLOBAL(selected_GLOBAL) + 0.032786.*temp_k_GLOBAL(selected_GLOBAL) - 14.8435 - 0.071692.*F1.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.0021487.*salinity_GLOBAL(selected_GLOBAL);
        K1(selected_GLOBAL)  = 10.^-pK1(selected_GLOBAL)...         % this is on the NBS scale
            ./fH(selected_GLOBAL);                    % convert to SWS scale (uncertain at low salinity_GLOBAL due to junction potential);
	    F2 = -129.24./temp_k_GLOBAL(selected_GLOBAL) + 1.4381;
	    pK2(selected_GLOBAL) = 2902.39./temp_k_GLOBAL(selected_GLOBAL) + 0.02379.*temp_k_GLOBAL(selected_GLOBAL) - 6.4980 - 0.3191.*F2.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.0198.*salinity_GLOBAL(selected_GLOBAL);
        K2(selected_GLOBAL)  = 10.^-pK2(selected_GLOBAL)...         % this is on the NBS scale
            ./fH(selected_GLOBAL);                    % convert to SWS scale (uncertain at low salinity_GLOBAL due to junction potential); 
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==10);
    if any(selected_GLOBAL)
        % From Lueker, Dickson, Keeling, 2000
	    % This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work. 
        % Mar. Chem. 70 (2000) 105-119
        % Total scale and kg-sw
        pK1(selected_GLOBAL) = 3633.86./temp_k_GLOBAL(selected_GLOBAL)-61.2172+9.6777.*log(temp_k_GLOBAL(selected_GLOBAL))-0.011555.*salinity_GLOBAL(selected_GLOBAL)+0.0001152.*salinity_GLOBAL(selected_GLOBAL).^2;
	    K1(selected_GLOBAL)  = 10.^-pK1(selected_GLOBAL)...           % this is on the total pH scale in mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                % convert to SWS pH scale
        pK2(selected_GLOBAL) = 471.78./temp_k_GLOBAL(selected_GLOBAL)+25.929 -3.16967.*log(temp_k_GLOBAL(selected_GLOBAL))-0.01781 .*salinity_GLOBAL(selected_GLOBAL)+0.0001122.*salinity_GLOBAL(selected_GLOBAL).^2;
	    K2(selected_GLOBAL)  = 10.^-pK2(selected_GLOBAL)...           % this is on the total pH scale in mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                % convert to SWS pH scale
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==11);
    if any(selected_GLOBAL)
	    % Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
	    % sigma for pK1 is reported to be 0.0056
	    % sigma for pK2 is reported to be 0.010
	    % This is from the abstract and pages 2536-2537
        pK1 =  -43.6977 - 0.0129037.*salinity_GLOBAL(selected_GLOBAL) + 1.364e-4.*salinity_GLOBAL(selected_GLOBAL).^2 + 2885.378./temp_k_GLOBAL(selected_GLOBAL) +  7.045159.*log(temp_k_GLOBAL(selected_GLOBAL));
        pK2 = -452.0940 + 13.142162.*salinity_GLOBAL(selected_GLOBAL) - 8.101e-4.*salinity_GLOBAL(selected_GLOBAL).^2 + 21263.61./temp_k_GLOBAL(selected_GLOBAL) + 68.483143.*log(temp_k_GLOBAL(selected_GLOBAL))...
				    + (-581.4428.*salinity_GLOBAL(selected_GLOBAL) + 0.259601.*salinity_GLOBAL(selected_GLOBAL).^2)./temp_k_GLOBAL(selected_GLOBAL) - 1.967035.*salinity_GLOBAL(selected_GLOBAL).*log(temp_k_GLOBAL(selected_GLOBAL));
	    K1(selected_GLOBAL) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	    K2(selected_GLOBAL) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==12);
    if any(selected_GLOBAL)
	    % Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
	    % Calculated from overdetermined WOCE-era field measurements 
	    % sigma for pK1 is reported to be 0.005
	    % sigma for pK2 is reported to be 0.008
	    % This is from page 1715
        pK1 =  6.359 - 0.00664.*salinity_GLOBAL(selected_GLOBAL) - 0.01322.*TempC(selected_GLOBAL) + 4.989e-5.*TempC(selected_GLOBAL).^2;
        pK2 =  9.867 - 0.01314.*salinity_GLOBAL(selected_GLOBAL) - 0.01904.*TempC(selected_GLOBAL) + 2.448e-5.*TempC(selected_GLOBAL).^2;
	    K1(selected_GLOBAL) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	    K2(selected_GLOBAL) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==13);
    if any(selected_GLOBAL)
        % From Millero 2006 work on pK1 and pK2 from titrations
	    % Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
        % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
	    pK1_0 = -126.34048 + 6320.813./temp_k_GLOBAL(selected_GLOBAL) + 19.568224*log(temp_k_GLOBAL(selected_GLOBAL));
	    A_1   = 13.4191*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.0331.*salinity_GLOBAL(selected_GLOBAL) - 5.33e-5.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B_1   = -530.123*salinity_GLOBAL(selected_GLOBAL).^0.5 - 6.103.*salinity_GLOBAL(selected_GLOBAL);
	    C_1   = -2.06950.*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK1(selected_GLOBAL)= A_1 + B_1./temp_k_GLOBAL(selected_GLOBAL) + C_1.*log(temp_k_GLOBAL(selected_GLOBAL)) + pK1_0; % pK1 sigma = 0.0054
        K1(selected_GLOBAL) = 10.^-(pK1(selected_GLOBAL));
	    pK2_0= -90.18333 + 5143.692./temp_k_GLOBAL(selected_GLOBAL) + 14.613358*log(temp_k_GLOBAL(selected_GLOBAL));	
	    A_2   = 21.0894*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.1248.*salinity_GLOBAL(selected_GLOBAL) - 3.687e-4.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B_2   = -772.483*salinity_GLOBAL(selected_GLOBAL).^0.5 - 20.051.*salinity_GLOBAL(selected_GLOBAL);
	    C_2   = -3.3336.*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK2(selected_GLOBAL)= A_2 + B_2./temp_k_GLOBAL(selected_GLOBAL) + C_2.*log(temp_k_GLOBAL(selected_GLOBAL)) + pK2_0; %pK2 sigma = 0.011
        K2(selected_GLOBAL) = 10.^-(pK2(selected_GLOBAL));
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==14);
    if any(selected_GLOBAL)
      % From Millero, 2010, also for estuarine use.
	    % Marine and Freshwater Research, v. 61, p. 139-142.
	    % Fits through compilation of real seawater titration results:
	    % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
	    % Constants for K's on the SWS;
	    % This is from page 141
	    pK10 = -126.34048 + 6320.813./temp_k_GLOBAL(selected_GLOBAL) + 19.568224.*log(temp_k_GLOBAL(selected_GLOBAL));
	    % This is from their table 2, page 140.
	    A1 = 13.4038.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.03206.*salinity_GLOBAL(selected_GLOBAL) - 5.242e-5.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B1 = -530.659.*salinity_GLOBAL(selected_GLOBAL).^0.5 - 5.8210.*salinity_GLOBAL(selected_GLOBAL);
	    C1 = -2.0664*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK1 = pK10 + A1 + B1./temp_k_GLOBAL(selected_GLOBAL) + C1.*log(temp_k_GLOBAL(selected_GLOBAL));
	    K1(selected_GLOBAL) = 10.^-pK1;
	    % This is from page 141
	    pK20 =  -90.18333 + 5143.692./temp_k_GLOBAL(selected_GLOBAL) + 14.613358.*log(temp_k_GLOBAL(selected_GLOBAL));
	    % This is from their table 3, page 140.
	    A2 = 21.3728.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.1218.*salinity_GLOBAL(selected_GLOBAL) - 3.688e-4.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B2 = -788.289.*salinity_GLOBAL(selected_GLOBAL).^0.5 - 19.189.*salinity_GLOBAL(selected_GLOBAL);
	    C2 = -3.374.*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK2 = pK20 + A2 + B2./temp_k_GLOBAL(selected_GLOBAL) + C2.*log(temp_k_GLOBAL(selected_GLOBAL));
	    K2(selected_GLOBAL) = 10.^-pK2;
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==15);
    % Added by J. C. Orr on 4 Dec 2016
    if any(selected_GLOBAL)
        % From Waters, Millero, and Woosley, 2014
	    % Mar. Chem., 165, 66-67, 2014
      % Corrigendum to "The free proton concentration scale for seawater pH".
	    % Effectively, this is an update of Millero (2010) formulation (which_k1_k2_constants_GLOBAL==14)
	    % Constants for K's on the SWS;
	    pK10 = -126.34048 + 6320.813./temp_k_GLOBAL(selected_GLOBAL) + 19.568224.*log(temp_k_GLOBAL(selected_GLOBAL));
	    A1 = 13.409160.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.031646.*salinity_GLOBAL(selected_GLOBAL) - 5.1895e-5.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B1 = -531.3642.*salinity_GLOBAL(selected_GLOBAL).^0.5 - 5.713.*salinity_GLOBAL(selected_GLOBAL);
	    C1 = -2.0669166.*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK1 = pK10 + A1 + B1./temp_k_GLOBAL(selected_GLOBAL) + C1.*log(temp_k_GLOBAL(selected_GLOBAL));
	    K1(selected_GLOBAL) = 10.^-pK1;
	    pK20 =  -90.18333 + 5143.692./temp_k_GLOBAL(selected_GLOBAL) + 14.613358.*log(temp_k_GLOBAL(selected_GLOBAL));
	    A2 = 21.225890.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.12450870.*salinity_GLOBAL(selected_GLOBAL) - 3.7243e-4.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B2 = -779.3444.*salinity_GLOBAL(selected_GLOBAL).^0.5 - 19.91739.*salinity_GLOBAL(selected_GLOBAL);
	    C2 = -3.3534679.*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK2 = pK20 + A2 + B2./temp_k_GLOBAL(selected_GLOBAL) + C2.*log(temp_k_GLOBAL(selected_GLOBAL));
	    K2(selected_GLOBAL) = 10.^-pK2;
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==16);
    % Added by J. D. Sharp on 9 Jul 2020
    if any(selected_GLOBAL)
        % From Sulpis et al, 2020
	    % Ocean Science Discussions, 16, 847-862
        % This study uses overdeterminations of the carbonate system to
        % iteratively fit K1 and K2
        pK1(selected_GLOBAL) = 8510.63./temp_k_GLOBAL(selected_GLOBAL)-172.4493+26.32996.*log(temp_k_GLOBAL(selected_GLOBAL))-0.011555.*salinity_GLOBAL(selected_GLOBAL)+0.0001152.*salinity_GLOBAL(selected_GLOBAL).^2;
	    K1(selected_GLOBAL)  = 10.^-pK1(selected_GLOBAL)...           % this is on the total pH scale in mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                % convert to SWS pH scale
        pK2(selected_GLOBAL) = 4226.23./temp_k_GLOBAL(selected_GLOBAL)-59.4636+9.60817.*log(temp_k_GLOBAL(selected_GLOBAL))-0.01781 .*salinity_GLOBAL(selected_GLOBAL)+0.0001122.*salinity_GLOBAL(selected_GLOBAL).^2;
	    K2(selected_GLOBAL)  = 10.^-pK2(selected_GLOBAL)...           % this is on the total pH scale in mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                % convert to SWS pH scale
    end
    
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==17);
    % Added by J. D. Sharp on 15 Feb 2021
    if any(selected_GLOBAL)
        % From Schockman & Byrne, 2021
	    % Geochimica et Cosmochimica Acta, in press
        % This study uses spectrophotometric pH measurements to determine
        % K1*K2 with unprecedented precision, and presents a new
        % parameterization for K2 based on these determinations
        % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
        pK10 = -126.34048 + 6320.813./temp_k_GLOBAL(selected_GLOBAL) + 19.568224.*log(temp_k_GLOBAL(selected_GLOBAL));
	    A1 = 13.568513.*salinity_GLOBAL(selected_GLOBAL).^0.5 + 0.031645.*salinity_GLOBAL(selected_GLOBAL) - 5.3834e-5.*salinity_GLOBAL(selected_GLOBAL).^2;
	    B1 = -539.2304.*salinity_GLOBAL(selected_GLOBAL).^0.5 - 5.635.*salinity_GLOBAL(selected_GLOBAL);
	    C1 = -2.0901396.*salinity_GLOBAL(selected_GLOBAL).^0.5;
	    pK1 = pK10 + A1 + B1./temp_k_GLOBAL(selected_GLOBAL) + C1.*log(temp_k_GLOBAL(selected_GLOBAL));
	    K1(selected_GLOBAL) = 10.^-pK1...               % this is on the total pH scale in mol/kg-sw
            ./SWStoTOT(selected_GLOBAL);                % convert to SWS pH scale
        % K2 is based on measurements of K1*K2:
        pK2 = 116.8067 - 3655.02./temp_k_GLOBAL(selected_GLOBAL) - 16.45817.*log(temp_k_GLOBAL(selected_GLOBAL)) + ...
            0.04523.*salinity_GLOBAL(selected_GLOBAL) - 0.615.*salinity_GLOBAL(selected_GLOBAL).^0.5 - 0.0002799.*salinity_GLOBAL(selected_GLOBAL).^2 + ...
            4.969.*(salinity_GLOBAL(selected_GLOBAL)./temp_k_GLOBAL(selected_GLOBAL));
        K2(selected_GLOBAL)  = 10.^-pK2...           % this is on the total pH scale in mol/kg-SW
            ./SWStoTOT(selected_GLOBAL);                % convert to SWS pH scale
    end
    
    %***************************************************************************
    %CorrectKsForPressureNow:
    % Currently: For which_k1_k2_constants_GLOBAL% = 1 to 7, all Ks (except KF and KS, which are on
    %       the free scale) are on the SWS scale.
    %       For which_k1_k2_constants_GLOBAL% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
    %       For which_k1_k2_constants_GLOBAL% = 8, K1, K2, and KW are on the "pH" pH scale
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
    %       Millero, selected_GLOBAL. J., Thermodynamics of the carbon dioxide system in the
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
    %       Millero, selected_GLOBAL. J., The thermodynamics of the carbon dioxide system
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
    % For which_k1_k2_constants_GLOBAL% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
    %       KP3, and KSi as for the other cases. Peng et al didn't consider the
    %       case of P different from 0. GEOSECS did consider pressure, but didn't
    %       include Phos, Si, or OH, so including the factors here won't matter.
    % For which_k1_k2_constants_GLOBAL% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
    %       and KW). The other aren't used (boron_concentration_GLOBAL = sulphate_concentration_GLOBAL = fluorine_concentration_GLOBAL = phosphate_GLOBAL = silicate_GLOBAL = 0.), so
    %       including the factors won't matter.
    %****************************************************************************
    %       deltaVs are in cm3/mole
    %       Kappas are in cm3/mole/bar
    %****************************************************************************
    
    %CorrectK1K2KBForPressure:
    deltaV    = nan(ntps,1); Kappa     = nan(ntps,1);
    lnK1fac   = nan(ntps,1); lnK2fac   = nan(ntps,1);
    lnKBfac   = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8);
    RR = (gas_constant.*temp_k_GLOBAL);
    if any(selected_GLOBAL)
        %***PressureEffectsOnK1inFreshWater:
        %               This is from Millero, 1983.
        deltaV(selected_GLOBAL)  = -30.54 + 0.1849 .*TempC(selected_GLOBAL) - 0.0023366.*TempC(selected_GLOBAL).^2;
        Kappa(selected_GLOBAL)   = (-6.22 + 0.1368 .*TempC(selected_GLOBAL) - 0.001233 .*TempC(selected_GLOBAL).^2)./1000;
        lnK1fac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        %***PressureEffectsOnK2inFreshWater:
        %               This is from Millero, 1983.
        deltaV(selected_GLOBAL)  = -29.81 + 0.115.*TempC(selected_GLOBAL) - 0.001816.*TempC(selected_GLOBAL).^2;
        Kappa(selected_GLOBAL)   = (-5.74 + 0.093.*TempC(selected_GLOBAL) - 0.001896.*TempC(selected_GLOBAL).^2)./1000;
        lnK2fac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        lnKBfac(selected_GLOBAL) = 0 ;%; this doesn't matter since boron_concentration_GLOBAL = 0 for this case
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
        %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
        %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
        %               but the fits are the same as those in
        %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
        %               who in turn quote Li, personal communication
        lnK1fac(selected_GLOBAL) = (24.2 - 0.085.*TempC(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        lnK2fac(selected_GLOBAL) = (16.4 - 0.04 .*TempC(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
        %               and matches the GEOSECS results
        lnKBfac(selected_GLOBAL) = (27.5 - 0.095.*TempC(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8);
    if any(selected_GLOBAL)
        %***PressureEffectsOnK1:
        %               These are from Millero, 1995.
        %               They are the same as Millero, 1979 and Millero, 1992.
        %               They are from data of Culberson and Pytkowicz, 1968.
        deltaV(selected_GLOBAL)  = -25.5 + 0.1271.*TempC(selected_GLOBAL);
        %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
        Kappa(selected_GLOBAL)   = (-3.08 + 0.0877.*TempC(selected_GLOBAL))./1000;
        %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
 	    lnK1fac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        %               The fits given in Millero, 1983 are somewhat different.
        
        %***PressureEffectsOnK2:
        %               These are from Millero, 1995.
        %               They are the same as Millero, 1979 and Millero, 1992.
        %               They are from data of Culberson and Pytkowicz, 1968.
        deltaV(selected_GLOBAL)  = -15.82 - 0.0219.*TempC(selected_GLOBAL);
        %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
        Kappa(selected_GLOBAL)   = (1.13 - 0.1475.*TempC(selected_GLOBAL))./1000;
        %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
	    lnK2fac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        %               The fit given in Millero, 1983 is different.
        %               Not by a lot for deltaV, but by much for Kappa. %
        
        %***PressureEffectsOnKB:
        %               This is from Millero, 1979.
        %               It is from data of Culberson and Pytkowicz, 1968.
        deltaV(selected_GLOBAL)  = -29.48 + 0.1622.*TempC(selected_GLOBAL) - 0.002608.*TempC(selected_GLOBAL).^2;
        %               Millero, 1983 has:
        %                 'deltaV = -28.56 + .1211.*temperature_in_GLOBAL - .000321.*temperature_in_GLOBAL.*temperature_in_GLOBAL
        %               Millero, 1992 has:
        %                 'deltaV = -29.48 + .1622.*temperature_in_GLOBAL + .295.*(Sali - 34.8)
        %               Millero, 1995 has:
        %                 'deltaV = -29.48 - .1622.*temperature_in_GLOBAL - .002608.*temperature_in_GLOBAL.*temperature_in_GLOBAL
        %                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979
        Kappa(selected_GLOBAL)   = -2.84./1000; % Millero, 1979
        %               Millero, 1992 and Millero, 1995 also have this.
        %                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
        %               Millero, 1983 has:
        %                 'Kappa = (-3 + .0427.*temperature_in_GLOBAL)./1000
        lnKBfac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
    end
    
    % CorrectKWForPressure:
    lnKWfac   = nan(ntps,1);
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8);
    if any(selected_GLOBAL)
        % PressureEffectsOnKWinFreshWater:
        %               This is from Millero, 1983.
        deltaV(selected_GLOBAL)  =  -25.6 + 0.2324.*TempC(selected_GLOBAL) - 0.0036246.*TempC(selected_GLOBAL).^2;
        Kappa(selected_GLOBAL)   = (-7.33 + 0.1368.*TempC(selected_GLOBAL) - 0.001233 .*TempC(selected_GLOBAL).^2)./1000;
 	    lnKWfac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
    
        %               NOTE the temperature dependence of KappaK1 and KappaKW
        %               for fresh water in Millero, 1983 are the same.
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=8);
    if any(selected_GLOBAL)
        % GEOSECS doesn't include OH term, so this won't matter.
        % Peng et al didn't include pressure, but here I assume that the KW correction
        %       is the same as for the other seawater cases.
        % PressureEffectsOnKW:
        %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
        deltaV(selected_GLOBAL)  = -20.02 + 0.1119.*TempC(selected_GLOBAL) - 0.001409.*TempC(selected_GLOBAL).^2;
        %               Millero, 1992 and Millero, 1995 have:
        Kappa(selected_GLOBAL)   = (-5.13 + 0.0794.*TempC(selected_GLOBAL))./1000; % Millero, 1983
        %               Millero, 1995 has this too, but Millero, 1992 is different.
	    lnKWfac(selected_GLOBAL) = (-deltaV(selected_GLOBAL) + 0.5.*Kappa(selected_GLOBAL).*Pbar(selected_GLOBAL)).*Pbar(selected_GLOBAL)./RR(selected_GLOBAL);
        %               Millero, 1979 does not list values for these.
    end
    
    % PressureEffectsOnKF:
    %       This is from Millero, 1995, which is the same as Millero, 1983.
    %       It is assumed that KF is on the free pH scale.
    deltaV = -9.78 - 0.009.*TempC - 0.000942.*TempC.^2;
    Kappa = (-3.91 + 0.054.*TempC)./1000;
    lnKFfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    % PressureEffectsOnKS:
    %       This is from Millero, 1995, which is the same as Millero, 1983.
    %       It is assumed that KS is on the free pH scale.
    deltaV = -18.03 + 0.0466.*TempC + 0.000316.*TempC.^2;
    Kappa = (-4.53 + 0.09.*TempC)./1000;
    lnKSfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    
    % CorrectKP1KP2KP3KSiForPressure:
    % These corrections don't matter for the GEOSECS choice (which_k1_k2_constants_GLOBAL% = 6) and
    %       the freshwater choice (which_k1_k2_constants_GLOBAL% = 8). For the Peng choice I assume
    %       that they are the same as for the other choices (which_k1_k2_constants_GLOBAL% = 1 to 5).
    % The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
    %       same as Millero, 1983.
    % PressureEffectsOnKP1:
    deltaV = -14.51 + 0.1211.*TempC - 0.000321.*TempC.^2;
    Kappa  = (-2.67 + 0.0427.*TempC)./1000;
    lnKP1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    % PressureEffectsOnKP2:
    deltaV = -23.12 + 0.1758.*TempC - 0.002647.*TempC.^2;
    Kappa  = (-5.15 + 0.09  .*TempC)./1000;
    lnKP2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    % PressureEffectsOnKP3:
    deltaV = -26.57 + 0.202 .*TempC - 0.003042.*TempC.^2;
    Kappa  = (-4.08 + 0.0714.*TempC)./1000;
    lnKP3fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    % PressureEffectsOnKSi:
    %  The only mention of this is Millero, 1995 where it is stated that the
    %    values have been estimated from the values of boric acid. HOWEVER,
    %    there is no listing of the values in the table.
    %    I used the values for boric acid from above.
    deltaV = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;
    Kappa  = -2.84./1000;
    lnKSifac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    % PressureEffectsOnKNH4: added by J. Sharp
    deltaV = -26.43 + 0.0889.*TempC - 0.000905.*TempC.^2;
    Kappa  = (-5.03 + 0.0814.*TempC)./1000;
    lnKNH4fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    % PressureEffectsOnKH2S: added by J. Sharp
    deltaV = -11.07 - 0.009.*TempC - 0.000942.*TempC.^2;
    Kappa  = (-2.89 + 0.054 .*TempC)./1000;
    lnKH2Sfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./(gas_constant.*temp_k_GLOBAL);
    
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
    SWStoTOT  = (1 + sulphate_concentration_GLOBAL./KS)./(1 + sulphate_concentration_GLOBAL./KS + fluorine_concentration_GLOBAL./KF);
    FREEtoTOT =  1 + sulphate_concentration_GLOBAL./KS;
    
    %  The values KS and KF are already pressure-corrected, so the pH scale
    %  conversions are now valid at pressure.
    
    % FindpHScaleConversionFactor:
    % this is the scale they will be put on
    pHfactor  = nan(ntps,1);
    selected_GLOBAL=(pH_scale==1); %Total
    pHfactor(selected_GLOBAL) = SWStoTOT(selected_GLOBAL);
    selected_GLOBAL=(pH_scale==2); %SWS, they are all on this now
    pHfactor(selected_GLOBAL) = 1;
    selected_GLOBAL=(pH_scale==3); %pHfree
    pHfactor(selected_GLOBAL) = SWStoTOT(selected_GLOBAL)./FREEtoTOT(selected_GLOBAL);
    selected_GLOBAL=(pH_scale==4); %pHNBS
    pHfactor(selected_GLOBAL) = fH(selected_GLOBAL);
    
    % ConvertFromSWSpHScaleToChosenScale:
    K1   = K1.* pHfactor;  K2   = K2.* pHfactor;
    KW   = KW.* pHfactor;  KB   = KB.* pHfactor;
    KP1  = KP1.*pHfactor;  KP2  = KP2.*pHfactor;
    KP3  = KP3.*pHfactor;  KSi  = KSi.*pHfactor;
    KNH4 = KNH4.*pHfactor; KH2S = KH2S.*pHfactor;
    



    Ks = containers.Map(["K0","K1","K2","KW","KB","KF","KS","KP1","KP2","KP3","KSi","KNH4","KH2S"], ...
                        {K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S});
end
