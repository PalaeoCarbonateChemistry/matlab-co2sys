classdef EquilibriumConstants
    methods (Static=true)
        %% Surface
        function k0 = calculate_surface_k0(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);
        
            % CalculateK0:
            % Weiss, R. selected., Marine Chemistry 2:203-215, 1974.
            TempK100  = temp_k./100;
            lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + salinity .*...
                (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
            K0   = exp(lnK0);                 % this is in mol/kg-SW/atm
            vCO2 = 32.3;                      % partial molal volume of CO2 (cm3 / mol)
                                              % from Weiss (1974, Appendix, paragraph 3)
            k0 = K0;
        end
        function [k1,k2] = calculate_surface_k1_and_k2(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);

            fH = calculate_fH(number_of_points,which_k1_k2,salinity,temp_k);
        
            % CalculateK1K2:
            logK1    = nan(number_of_points,1); lnK1     = nan(number_of_points,1);
            pK1      = nan(number_of_points,1); K1       = nan(number_of_points,1);
            logK2    = nan(number_of_points,1); lnK2     = nan(number_of_points,1);
            pK2      = nan(number_of_points,1); K2       = nan(number_of_points,1);

            
            selected=(which_k1_k2==1);
            if any(selected)
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
                lnK1(selected) = 2.83655 - 2307.1266./temp_k(selected) - 1.5529413.*log_temp_k(selected) +...
                    (-0.20760841 - 4.0484./temp_k(selected)).*sqrt(salinity(selected)) + 0.08468345.*salinity(selected) -...
                    0.00654208.*sqrt(salinity(selected)).*salinity(selected);
                K1(selected) = exp(lnK1(selected))...            % this is on the total pH scale in mol/kg-H2O
                    .*(1 - 0.001005.*salinity(selected))...    % convert to mol/kg-SW
                    ./SWStoTOT(selected);                 % convert to SWS pH scale
                % This is eq. 30 on p. 254 and what they use in their abstract:
                lnK2(selected) = -9.226508 - 3351.6106./temp_k(selected) - 0.2005743.*log_temp_k(selected) +...
                    (-0.106901773 - 23.9722./temp_k(selected)).*sqrt(salinity(selected)) + 0.1130822.*salinity(selected) -...
                    0.00846934.*sqrt(salinity(selected)).*salinity(selected);
                K2(selected) = exp(lnK2(selected))...            % this is on the total pH scale in mol/kg-H2O
                    .*(1 - 0.001005.*salinity(selected))...    % convert to mol/kg-SW
                    ./SWStoTOT(selected);                 % convert to SWS pH scale
            end
            selected=(which_k1_k2==2);
            if any(selected)
                % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
                % The 2s precision in pK1 is .011, or 2.5% in K1.
                % The 2s precision in pK2 is .02, or 4.5% in K2.
                % This is in Table 5 on p. 1652 and what they use in the abstract:
                pK1(selected) = 812.27./temp_k(selected) + 3.356 - 0.00171.*salinity(selected).*log_temp_k(selected)...
                    + 0.000091.*salinity(selected).^2;
                K1(selected) = 10.^(-pK1(selected)); % this is on the SWS pH scale in mol/kg-SW
                % This is in Table 5 on p. 1652 and what they use in the abstract:
                pK2(selected) = 1450.87./temp_k(selected) + 4.604 - 0.00385.*salinity(selected).*log_temp_k(selected)...
                    + 0.000182.*salinity(selected).^2;
                K2(selected) = 10.^(-pK2(selected)); % this is on the SWS pH scale in mol/kg-SW
            end
            selected=(which_k1_k2==3);
            if any(selected)
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
                pK1(selected) = 851.4./temp_k(selected) + 3.237 - 0.0106.*salinity(selected) + 0.000105.*salinity(selected).^2;
                K1(selected) = 10.^(-pK1(selected)); % this is on the SWS pH scale in mol/kg-SW
                % This is from Table 4 on p. 1739.
                pK2(selected) = -3885.4./temp_k(selected) + 125.844 - 18.141.*log_temp_k(selected)...
                    - 0.0192.*salinity(selected) + 0.000132.*salinity(selected).^2;
                K2(selected) = 10.^(-pK2(selected)); % this is on the SWS pH scale in mol/kg-SW
            end
            selected=(which_k1_k2==4);
            if any(selected)
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
                pK1(selected) = 3670.7./temp_k(selected) - 62.008 + 9.7944.*log_temp_k(selected)...
                         - 0.0118.*salinity(selected) + 0.000116.*salinity(selected).^2;
                K1(selected) = 10.^(-pK1(selected)); % this is on the SWS pH scale in mol/kg-SW
                % This is in Table 4 on p. 1739.
                pK2(selected) = 1394.7./temp_k(selected) + 4.777 - 0.0184.*salinity(selected) + 0.000118.*salinity(selected).^2;
                K2(selected) = 10.^(-pK2(selected)); % this is on the SWS pH scale in mol/kg-SW
            end
            selected=(which_k1_k2==5);
            if any(selected)
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
                pK1(selected) = 845./temp_k(selected) + 3.248 - 0.0098.*salinity(selected) + 0.000087.*salinity(selected).^2;
                K1(selected) = 10.^(-pK1(selected)); % this is on the SWS pH scale in mol/kg-SW
                % This is in Table 5 on p. 1740.
                pK2(selected) = 1377.3./temp_k(selected) + 4.824 - 0.0185.*salinity(selected) + 0.000122.*salinity(selected).^2;
                K2(selected) = 10.^(-pK2(selected)); % this is on the SWS pH scale in mol/kg-SW
            end
            selected=(which_k1_k2==6 | which_k1_k2==7);
            if any(selected)
                % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
                % Limnology and Oceanography, 18(6):897-907, 1973.
	            % I.e., these are the original Mehrbach dissociation constants.
                % The 2s precision in pK1 is .005, or 1.2% in K1.
                % The 2s precision in pK2 is .008, or 2% in K2.
                pK1(selected) = - 13.7201 + 0.031334.*temp_k(selected) + 3235.76./temp_k(selected)...
                    + 1.3e-5*salinity(selected).*temp_k(selected) - 0.1032.*salinity(selected).^0.5;
                K1(selected) = 10.^(-pK1(selected))...         % this is on the NBS scale
                    ./fH(selected);                     % convert to SWS scale
                pK2(selected) = 5371.9645 + 1.671221.*temp_k(selected) + 0.22913.*salinity(selected) + 18.3802.*log10(salinity(selected))...
                         - 128375.28./temp_k(selected) - 2194.3055.*log10(temp_k(selected)) - 8.0944e-4.*salinity(selected).*temp_k(selected)...
                         - 5617.11.*log10(salinity(selected))./temp_k(selected) + 2.136.*salinity(selected)./temp_k(selected); % pK2 is not defined for salinity=0, since log10(0)=-inf
                K2(selected) = 10.^(-pK2(selected))...         % this is on the NBS scale
                    ./fH(selected);                     % convert to SWS scale
            end
            selected=(which_k1_k2==8);
            if any(selected)	
	            % PURE WATER CASE
                % Millero, selected. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
                % K1 from refit data from Harned and Davis,
                % J American Chemical Society, 65:2030-2037, 1943.
                % K2 from refit data from Harned and Scholes,
                % J American Chemical Society, 43:1706-1709, 1941.
	            % This is only to be used for salinity=0 water (note the absence of S in the below formulations)
                % These are the thermodynamic Constants:
                lnK1(selected) = 290.9097 - 14554.21./temp_k(selected) - 45.0575.*log_temp_k(selected);
                K1(selected) = exp(lnK1(selected));
                lnK2(selected) = 207.6548 - 11843.79./temp_k(selected) - 33.6485.*log_temp_k(selected);
                K2(selected) = exp(lnK2(selected));
            end
            selected=(which_k1_k2==9);
            if any(selected)
                % From Cai and Wang 1998, for estuarine use.
	            % Data used in this work is from:
	            % K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	            % K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	            % Sigma of residuals between fits and above data: ±0.015, +0.040 for K1 and K2, respectively.
	            % salinity 0-40, Temp 0.2-30
              % Limnol. Oceanogr. 43(4) (1998) 657-668
	            % On the NBS scale
	            % Their check values for F1 don't work out, not sure if this was correctly published...
	            F1 = 200.1./temp_k(selected) + 0.3220;
	            pK1(selected) = 3404.71./temp_k(selected) + 0.032786.*temp_k(selected) - 14.8435 - 0.071692.*F1.*salinity(selected).^0.5 + 0.0021487.*salinity(selected);
                K1(selected)  = 10.^-pK1(selected)...         % this is on the NBS scale
                    ./fH(selected);                    % convert to SWS scale (uncertain at low salinity due to junction potential);
	            F2 = -129.24./temp_k(selected) + 1.4381;
	            pK2(selected) = 2902.39./temp_k(selected) + 0.02379.*temp_k(selected) - 6.4980 - 0.3191.*F2.*salinity(selected).^0.5 + 0.0198.*salinity(selected);
                K2(selected)  = 10.^-pK2(selected)...         % this is on the NBS scale
                    ./fH(selected);                    % convert to SWS scale (uncertain at low salinity due to junction potential); 
            end
            selected=(which_k1_k2==10);
            if any(selected)
                % From Lueker, Dickson, Keeling, 2000
	            % This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work. 
                % Mar. Chem. 70 (2000) 105-119
                % Total scale and kg-sw
                pK1(selected) = 3633.86./temp_k(selected)-61.2172+9.6777.*log(temp_k(selected))-0.011555.*salinity(selected)+0.0001152.*salinity(selected).^2;
	            K1(selected)  = 10.^-pK1(selected)...           % this is on the total pH scale in mol/kg-SW
                    ./SWStoTOT(selected);                % convert to SWS pH scale
                pK2(selected) = 471.78./temp_k(selected)+25.929 -3.16967.*log(temp_k(selected))-0.01781 .*salinity(selected)+0.0001122.*salinity(selected).^2;
	            K2(selected)  = 10.^-pK2(selected)...           % this is on the total pH scale in mol/kg-SW
                    ./SWStoTOT(selected);                % convert to SWS pH scale
            end
            selected=(which_k1_k2==11);
            if any(selected)
	            % Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
	            % sigma for pK1 is reported to be 0.0056
	            % sigma for pK2 is reported to be 0.010
	            % This is from the abstract and pages 2536-2537
                pK1 =  -43.6977 - 0.0129037.*salinity(selected) + 1.364e-4.*salinity(selected).^2 + 2885.378./temp_k(selected) +  7.045159.*log(temp_k(selected));
                pK2 = -452.0940 + 13.142162.*salinity(selected) - 8.101e-4.*salinity(selected).^2 + 21263.61./temp_k(selected) + 68.483143.*log(temp_k(selected))...
				            + (-581.4428.*salinity(selected) + 0.259601.*salinity(selected).^2)./temp_k(selected) - 1.967035.*salinity(selected).*log(temp_k(selected));
	            K1(selected) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	            K2(selected) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
            end
            selected=(which_k1_k2==12);
            if any(selected)
	            % Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
	            % Calculated from overdetermined WOCE-era field measurements 
	            % sigma for pK1 is reported to be 0.005
	            % sigma for pK2 is reported to be 0.008
	            % This is from page 1715
                pK1 =  6.359 - 0.00664.*salinity(selected) - 0.01322.*temp_c(selected) + 4.989e-5.*temp_c(selected).^2;
                pK2 =  9.867 - 0.01314.*salinity(selected) - 0.01904.*temp_c(selected) + 2.448e-5.*temp_c(selected).^2;
	            K1(selected) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	            K2(selected) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
            end
            selected=(which_k1_k2==13);
            if any(selected)
                % From Millero 2006 work on pK1 and pK2 from titrations
	            % Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
                % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
	            pK1_0 = -126.34048 + 6320.813./temp_k(selected) + 19.568224*log(temp_k(selected));
	            A_1   = 13.4191*salinity(selected).^0.5 + 0.0331.*salinity(selected) - 5.33e-5.*salinity(selected).^2;
	            B_1   = -530.123*salinity(selected).^0.5 - 6.103.*salinity(selected);
	            C_1   = -2.06950.*salinity(selected).^0.5;
	            pK1(selected)= A_1 + B_1./temp_k(selected) + C_1.*log(temp_k(selected)) + pK1_0; % pK1 sigma = 0.0054
                K1(selected) = 10.^-(pK1(selected));
	            pK2_0= -90.18333 + 5143.692./temp_k(selected) + 14.613358*log(temp_k(selected));	
	            A_2   = 21.0894*salinity(selected).^0.5 + 0.1248.*salinity(selected) - 3.687e-4.*salinity(selected).^2;
	            B_2   = -772.483*salinity(selected).^0.5 - 20.051.*salinity(selected);
	            C_2   = -3.3336.*salinity(selected).^0.5;
	            pK2(selected)= A_2 + B_2./temp_k(selected) + C_2.*log(temp_k(selected)) + pK2_0; %pK2 sigma = 0.011
                K2(selected) = 10.^-(pK2(selected));
            end
            selected=(which_k1_k2==14);
            if any(selected)
              % From Millero, 2010, also for estuarine use.
	            % Marine and Freshwater Research, v. 61, p. 139-142.
	            % Fits through compilation of real seawater titration results:
	            % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
	            % Constants for K's on the SWS;
	            % This is from page 141
	            pK10 = -126.34048 + 6320.813./temp_k(selected) + 19.568224.*log(temp_k(selected));
	            % This is from their table 2, page 140.
	            A1 = 13.4038.*salinity(selected).^0.5 + 0.03206.*salinity(selected) - 5.242e-5.*salinity(selected).^2;
	            B1 = -530.659.*salinity(selected).^0.5 - 5.8210.*salinity(selected);
	            C1 = -2.0664*salinity(selected).^0.5;
	            pK1 = pK10 + A1 + B1./temp_k(selected) + C1.*log(temp_k(selected));
	            K1(selected) = 10.^-pK1;
	            % This is from page 141
	            pK20 =  -90.18333 + 5143.692./temp_k(selected) + 14.613358.*log(temp_k(selected));
	            % This is from their table 3, page 140.
	            A2 = 21.3728.*salinity(selected).^0.5 + 0.1218.*salinity(selected) - 3.688e-4.*salinity(selected).^2;
	            B2 = -788.289.*salinity(selected).^0.5 - 19.189.*salinity(selected);
	            C2 = -3.374.*salinity(selected).^0.5;
	            pK2 = pK20 + A2 + B2./temp_k(selected) + C2.*log(temp_k(selected));
	            K2(selected) = 10.^-pK2;
            end
            selected=(which_k1_k2==15);
            % Added by J. C. Orr on 4 Dec 2016
            if any(selected)
                % From Waters, Millero, and Woosley, 2014
	            % Mar. Chem., 165, 66-67, 2014
              % Corrigendum to "The free proton concentration scale for seawater pH".
	            % Effectively, this is an update of Millero (2010) formulation (which_k1_k2==14)
	            % Constants for K's on the SWS;
	            pK10 = -126.34048 + 6320.813./temp_k(selected) + 19.568224.*log(temp_k(selected));
	            A1 = 13.409160.*salinity(selected).^0.5 + 0.031646.*salinity(selected) - 5.1895e-5.*salinity(selected).^2;
	            B1 = -531.3642.*salinity(selected).^0.5 - 5.713.*salinity(selected);
	            C1 = -2.0669166.*salinity(selected).^0.5;
	            pK1 = pK10 + A1 + B1./temp_k(selected) + C1.*log(temp_k(selected));
	            K1(selected) = 10.^-pK1;
	            pK20 =  -90.18333 + 5143.692./temp_k(selected) + 14.613358.*log(temp_k(selected));
	            A2 = 21.225890.*salinity(selected).^0.5 + 0.12450870.*salinity(selected) - 3.7243e-4.*salinity(selected).^2;
	            B2 = -779.3444.*salinity(selected).^0.5 - 19.91739.*salinity(selected);
	            C2 = -3.3534679.*salinity(selected).^0.5;
	            pK2 = pK20 + A2 + B2./temp_k(selected) + C2.*log(temp_k(selected));
	            K2(selected) = 10.^-pK2;
            end
            selected=(which_k1_k2==16);
            % Added by J. D. Sharp on 9 Jul 2020
            if any(selected)
                % From Sulpis et al, 2020
	            % Ocean Science Discussions, 16, 847-862
                % This study uses overdeterminations of the carbonate system to
                % iteratively fit K1 and K2
                pK1(selected) = 8510.63./temp_k(selected)-172.4493+26.32996.*log(temp_k(selected))-0.011555.*salinity(selected)+0.0001152.*salinity(selected).^2;
	            K1(selected)  = 10.^-pK1(selected)...           % this is on the total pH scale in mol/kg-SW
                    ./SWStoTOT(selected);                % convert to SWS pH scale
                pK2(selected) = 4226.23./temp_k(selected)-59.4636+9.60817.*log(temp_k(selected))-0.01781 .*salinity(selected)+0.0001122.*salinity(selected).^2;
	            K2(selected)  = 10.^-pK2(selected)...           % this is on the total pH scale in mol/kg-SW
                    ./SWStoTOT(selected);                % convert to SWS pH scale
            end
            
            selected=(which_k1_k2==17);
            % Added by J. D. Sharp on 15 Feb 2021
            if any(selected)
                % From Schockman & Byrne, 2021
	            % Geochimica et Cosmochimica Acta, in press
                % This study uses spectrophotometric pH measurements to determine
                % K1*K2 with unprecedented precision, and presents a new
                % parameterization for K2 based on these determinations
                % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
                pK10 = -126.34048 + 6320.813./temp_k(selected) + 19.568224.*log(temp_k(selected));
	            A1 = 13.568513.*salinity(selected).^0.5 + 0.031645.*salinity(selected) - 5.3834e-5.*salinity(selected).^2;
	            B1 = -539.2304.*salinity(selected).^0.5 - 5.635.*salinity(selected);
	            C1 = -2.0901396.*salinity(selected).^0.5;
	            pK1 = pK10 + A1 + B1./temp_k(selected) + C1.*log(temp_k(selected));
	            K1(selected) = 10.^-pK1...               % this is on the total pH scale in mol/kg-sw
                    ./SWStoTOT(selected);                % convert to SWS pH scale
                % K2 is based on measurements of K1*K2:
                pK2 = 116.8067 - 3655.02./temp_k(selected) - 16.45817.*log(temp_k(selected)) + ...
                    0.04523.*salinity(selected) - 0.615.*salinity(selected).^0.5 - 0.0002799.*salinity(selected).^2 + ...
                    4.969.*(salinity(selected)./temp_k(selected));
                K2(selected)  = 10.^-pK2...           % this is on the total pH scale in mol/kg-SW
                    ./SWStoTOT(selected);                % convert to SWS pH scale
            end
            k1 = K1;
            k2 = K2;
        end
        function kw = calculate_surface_kw(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);
            
            lnKW = nan(number_of_points,1); KW = nan(number_of_points,1);
            selected=(which_k1_k2==7);
            if any(selected)
                % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
                lnKW(selected) = 148.9802 - 13847.26./temp_k(selected) - 23.6521.*log_temp_k(selected) +...
                    (-79.2447 + 3298.72./temp_k(selected) + 12.0408.*log_temp_k(selected)).*...
                    sqrt(salinity(selected)) - 0.019813.*salinity(selected);
            end
            selected=(which_k1_k2==8);
            if any(selected)
                % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
                % refit data of Harned and Owen, The Physical Chemistry of
                % Electrolyte Solutions, 1958
                lnKW(selected) = 148.9802 - 13847.26./temp_k(selected) - 23.6521.*log_temp_k(selected);
            end
            selected=(which_k1_k2~=6 & which_k1_k2~=7 & which_k1_k2~=8);
            if any(selected)
                % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
                % his check value of 1.6 umol/kg-SW should be 6.2
                lnKW(selected) = 148.9802 - 13847.26./temp_k(selected) - 23.6521.*log_temp_k(selected) +...
                    (-5.977 + 118.67./temp_k(selected) + 1.0495.*log_temp_k(selected)).*...
                    sqrt(salinity(selected)) - 0.01615.*salinity(selected);
            end
            KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
            selected=(which_k1_k2==6);
            if any(selected)
                KW(selected) = 0; % GEOSECS doesn't include OH effects
            end
            kw = KW;
        end
        function kb = calculate_surface_kb(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);
            fH = calculate_fH(number_of_points,which_k1_k2,salinity,temp_k);
            
            % CalculateKB:
            KB      = nan(number_of_points,1); logKB   = nan(number_of_points,1);
            lnKBtop = nan(number_of_points,1); lnKB    = nan(number_of_points,1);
            selected=(which_k1_k2==8); % Pure water case
            if any(selected)
                KB(selected) = 0;
            end
            selected=(which_k1_k2==6 | which_k1_k2==7);
            if any(selected)
                % This is for GEOSECS and Peng et al.
                % Lyman, John, UCLA Thesis, 1957
                % fit by Li et al, JGR 74:5507-5525, 1969:
                logKB(selected) = -9.26 + 0.00886.*salinity(selected) + 0.01.*temp_c(selected);
                KB(selected) = 10.^(logKB(selected))...  % this is on the NBS scale
                    ./fH(selected);               % convert to the SWS scale
            end
            selected=(which_k1_k2~=6 & which_k1_k2~=7 & which_k1_k2~=8);
            if any(selected)
                % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
                lnKBtop(selected) = -8966.9 - 2890.53.*sqrt(salinity(selected)) - 77.942.*salinity(selected) +...
                    1.728.*sqrt(salinity(selected)).*salinity(selected) - 0.0996.*salinity(selected).^2;
                lnKB(selected) = lnKBtop(selected)./temp_k(selected) + 148.0248 + 137.1942.*sqrt(salinity(selected)) +...
                    1.62142.*salinity(selected) + (-24.4344 - 25.085.*sqrt(salinity(selected)) - 0.2474.*...
                    salinity(selected)).*log_temp_k(selected) + 0.053105.*sqrt(salinity(selected)).*temp_k(selected);
                KB(selected) = exp(lnKB(selected))...    % this is on the total pH scale in mol/kg-SW
                    ./SWStoTOT(selected);         % convert to SWS pH scale
            end
            kb = KB;    
        end
        function [kp1,kp2,kp3,ksi] = calculate_surface_kp1_and_kp2_and_kp3_and_ksi(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);
            IonS         = 19.924 .* salinity ./ (1000 - 1.005   .* salinity);
            fH = calculate_fH(number_of_points,which_k1_k2,salinity,temp_k);
            
            % CalculateKP1KP2KP3KSi:
            KP1      = nan(number_of_points,1); KP2      = nan(number_of_points,1);
            KP3      = nan(number_of_points,1); KSi      = nan(number_of_points,1);
            lnKP1    = nan(number_of_points,1); lnKP2    = nan(number_of_points,1);
            lnKP3    = nan(number_of_points,1); lnKSi    = nan(number_of_points,1);
            selected=(which_k1_k2==7);
            if any(selected)
                KP1(selected) = 0.02;
                % Peng et al don't include the contribution from this term,
                % but it is so small it doesn't contribute. It needs to be
                % kept so that the routines work ok.
                % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
                % Limnology and Oceanography 12:243-252, 1967:
                % these are only for sals 33 to 36 and are on the NBS scale
                KP2(selected) = exp(-9.039 - 1450./temp_k(selected))... % this is on the NBS scale
                    ./fH(selected);                          % convert to SWS scale
                KP3(selected) = exp(4.466 - 7276./temp_k(selected))...  % this is on the NBS scale
                    ./fH(selected);                          % convert to SWS scale
                % Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
                % The Chemical Society (London), Special Publ. 17:751, 1964:
                KSi(selected) = 0.0000000004...              % this is on the NBS scale
                    ./fH(selected);                          % convert to SWS scale
                
            end
            selected=(which_k1_k2==6 | which_k1_k2==8);
            if any(selected)
                KP1(selected) = 0; KP2(selected) = 0; KP3(selected) = 0; KSi(selected) = 0;
                % Neither the GEOSECS choice nor the freshwater choice
                % include contributions from phosphate or silicate.
            end
            selected=(which_k1_k2~=6 & which_k1_k2~=7 & which_k1_k2~=8);
            if any(selected)
                % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
                % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
                % KSi was given on the SWS pH scale in molal units.
                lnKP1(selected) = -4576.752./temp_k(selected) + 115.54 - 18.453.*log_temp_k(selected) + (-106.736./temp_k(selected) +...
                    0.69171).*sqrt(salinity(selected)) + (-0.65643./temp_k(selected) - 0.01844).*salinity(selected);
                KP1(selected) = exp(lnKP1(selected));
                lnKP2(selected) = -8814.715./temp_k(selected) + 172.1033 - 27.927.*log_temp_k(selected) + (-160.34./temp_k(selected) +...
                    1.3566).*sqrt(salinity(selected)) + (0.37335./temp_k(selected) - 0.05778).*salinity(selected);
                KP2(selected) = exp(lnKP2(selected));
                lnKP3(selected) = -3070.75./temp_k(selected) - 18.126 + (17.27039./temp_k(selected) + 2.81197).*sqrt(salinity(selected)) +...
                    (-44.99486./temp_k(selected) - 0.09984).*salinity(selected);
                KP3(selected) = exp(lnKP3(selected));
                lnKSi(selected) = -8904.2./temp_k(selected) + 117.4 - 19.334.*log_temp_k(selected) + (-458.79./temp_k(selected) +...
                    3.5913).*sqrt(IonS(selected)) + (188.74./temp_k(selected) - 1.5998).*IonS(selected) +...
                    (-12.1652./temp_k(selected) + 0.07871).*IonS(selected).^2;
                KSi(selected) = exp(lnKSi(selected))...                % this is on the SWS pH scale in mol/kg-H2O
                    .*(1 - 0.001005.*salinity(selected));        % convert to mol/kg-SW
            end

            kp1 = KP1;
            kp2 = KP2;
            kp3 = KP3;
            ksi = KSi;
        end
        function [knh4,kh2s] = calculate_surface_knh4_and_kh2s(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);
            
            % Calculate KNH4 and KH2S: added by J. Sharp
            KNH4           = nan(number_of_points,1); KH2S       = nan(number_of_points,1);
            PKNH4expCW     = nan(number_of_points,1); lnKH2S     = nan(number_of_points,1);
            selected=(which_k1_k2==6 | which_k1_k2==7 | which_k1_k2==8); % GEOSECS or freshwater cases
            if any(selected)
                KNH4(selected) = 0;
                KH2S(selected) = 0;
            end
            selected=(which_k1_k2~=6 & which_k1_k2~=7 & which_k1_k2~=8); % All other cases
            if any(selected)
            % Ammonia dissociation constant from Yao and Millero (1995)
            %   KNH4(selected) = (exp(-6285.33./temp_k(selected)+0.0001635.*temp_k(selected)-0.25444+...
            %             (0.46532-123.7184./temp_k(selected)).*salinity(selected).^0.5+(-0.01992+...
            %             3.17556./temp_k(selected)).*salinity(selected)))...
            %             ./SWStoTOT(selected);                    % convert to SWS pH scale
            % Ammonia dissociation constant from Clegg and Whitfield (1995)
              PKNH4(selected) = 9.244605-2729.33.*(1./298.15-1./temp_k(selected)) +...
                      (0.04203362-11.24742./temp_k(selected)).*salinity(selected).^0.25+...  % added missing (selected) index on salinity // MPH
                      (-13.6416+1.176949.*temp_k(selected).^0.5-...
                      0.02860785.*temp_k(selected)+545.4834./temp_k(selected)).*salinity(selected).^0.5+...
                      (-0.1462507+0.0090226468.*temp_k(selected).^0.5-...
                      0.0001471361.*temp_k(selected)+10.5425./temp_k(selected)).*salinity(selected).^1.5+...
                      (0.004669309-0.0001691742.*temp_k(selected).^0.5-...
                      0.5677934./temp_k(selected)).*salinity(selected).^2+...
                      (-2.354039E-05+0.009698623./temp_k(selected)).*salinity(selected).^2.5;
              KNH4(selected)  = 10.^-PKNH4(selected);                    % total scale, mol/kg-H2O
              KNH4(selected)  = KNH4(selected).*(1-0.001005.*salinity(selected)); % mol/kg-SW
              KNH4(selected)  = KNH4(selected)./SWStoTOT(selected);             % converts to SWS pH scale
            
            % First hydrogen sulfide dissociation constant from Millero et al. (1988)
              KH2S(selected)  = (exp(225.838-13275.3./temp_k(selected)-34.6435.*log(temp_k(selected))+...
                          0.3449.*salinity(selected).^0.5-0.0274.*salinity(selected)))...
                          ./SWStoTOT(selected);                    % convert to SWS pH scale
            
            end
            knh4 = KNH4;
            kh2s = KH2S;
        end
        
        function ks = calculate_surface_ks(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);

            IonS         = 19.924 .* salinity ./ (1000 - 1.005   .* salinity);

            % CalculateKS:
            lnKS   = nan(number_of_points,1); pKS  = nan(number_of_points,1); KS   = nan(number_of_points,1);
            logKS0 = nan(number_of_points,1); logKSK0 = nan(number_of_points,1);
            selected=(which_kso4==1);
            if any(selected)
                % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
                % The goodness of fit is .021.
                % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
                % TYPO on p. 121: the constant e9 should be e8.
                % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
              lnKS(selected) = -4276.1./temp_k(selected) + 141.328 - 23.093.*log_temp_k(selected) +...             
                  (-13856./temp_k(selected) + 324.57 - 47.986.*log_temp_k(selected)).*sqrt(IonS(selected)) +...     
                  (35474./temp_k(selected) - 771.54 + 114.723.*log_temp_k(selected)).*IonS(selected) +...           
                  (-2698./temp_k(selected)).*sqrt(IonS(selected)).*IonS(selected) + (1776./temp_k(selected)).*IonS(selected).^2; 
	            KS(selected) = exp(lnKS(selected))...            % this is on the free pH scale in mol/kg-H2O
                    .* (1 - 0.001005 .* salinity(selected));   % convert to mol/kg-SW
            end
            selected=(which_kso4==2);
            if any(selected)
                % Khoo et al, Analytical Chemistry, 49(1):29-34, 1977
                % KS was found by titrations with a hydrogen electrode
                % of artificial seawater containing sulfate (but without selected)
                % at 3 salinities from 20 to 45 and artificial seawater NOT
                % containing sulfate (nor selected) at 16 salinities from 15 to 45,
                % both at temperatures from 5 to 40 deg C.
                % KS is on the Free pH scale (inherently so).
                % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
                % He finds log(beta) which = my pKS;
                % his beta is an association constant.
                % The rms error is .0021 in pKS, or about .5% in KS.
                % This is equation 20 on p. 33:
                pKS(selected) = 647.59 ./ temp_k(selected) - 6.3451 + 0.019085.*temp_k(selected) - 0.5208.*sqrt(IonS(selected));
                KS(selected) = 10.^(-pKS(selected))...          % this is on the free pH scale in mol/kg-H2O
                    .* (1 - 0.001005.*salinity(selected));    % convert to mol/kg-SW
            end
            selected=(which_kso4==3);
            if any(selected)
                % Waters and Millero, Marine Chemistry, 149: 8-22, 2013, with corrections from
                % Waters et al, Marine Chemistry, 165: 66-67, 2014
                logKS0(selected) = 562.69486 - 102.5154.*log_temp_k(selected) - 0.0001117033.*temp_k(selected).*temp_k(selected) + ...
                    0.2477538.*temp_k(selected) - 13273.76./temp_k(selected);
                logKSK0(selected) = (4.24666 - 0.152671.*temp_k(selected) + 0.0267059.*temp_k(selected).*log_temp_k(selected) - 0.000042128.*temp_k(selected).*temp_k(selected)).*salinity(selected).^0.5 + ...
                    (0.2542181 - 0.00509534.*temp_k(selected) + 0.00071589.*temp_k(selected).*log_temp_k(selected)).*salinity(selected) + (-0.00291179 + 0.0000209968.*temp_k(selected)).*salinity(selected).^1.5 + ...
                    -0.0000403724.*salinity(selected).^2;
                KS(selected) = ((10.^(logKSK0(selected))).*(10.^logKS0(selected))) ... % this is on the free pH scale in mol/kg-H2O
                    .* (1 - 0.001005.*salinity(selected));                    % convert to mol/kg-SW
            end
            ks = KS;
        end
        function kf = calculate_surface_kf(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2)
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);

            IonS         = 19.924 .* salinity ./ (1000 - 1.005   .* salinity);
            
            % CalculateKF:
            KF = NaN(number_of_points, 1);  % added preallocation here and selected-indexing below // MPH
            selected=(which_kf==1);
            if any(selected)
                % Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
                lnKF = 1590.2./temp_k - 12.641 + 1.525.*IonS.^0.5;
                KF(selected)   = exp(lnKF(selected))...                 % this is on the free pH scale in mol/kg-H2O
                    .*(1 - 0.001005.*salinity(selected));          % convert to mol/kg-SW
            end
            selected=(which_kf==2);
            if any(selected)
                % Perez and Fraga 1987 (to be used for S: 10-40, T: 9-33)
                % P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
                lnKF = 874./temp_k - 9.68 + 0.111.*salinity.^0.5;
                KF(selected)   = exp(lnKF(selected));                   % this is on the free pH scale in mol/kg-SW
            end
            kf = KF;
        end
        
        %% Pressure Corrections
        function k0_pressure_correction = calculate_pressure_correction_k0(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar,co2_correction)
            if co2_correction == 0
                co2_pressure_correction = 1; % Set pressure correction to 1
            elseif co2_correction == 1
                co2_pressure_correction = exp((-pressure_bar).*vCO2./(gas_constant.*temp_k)); % Calculate pressure correction to K0
            else
                disp('co2_press must be set to either 0 or 1'); % Display error message
            end         
            k0_pressure_correction = co2_pressure_correction; % this is in mol/kg-SW/atm
        end
        function [k1_pressure_correction,...
                  k2_pressure_correction,...
                  kb_pressure_correction] = calculate_pressure_correction_k1_and_k2_and_kb(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            
            temp_k    = temp_c + 273.15;
            log_temp_k = log(temp_k);
        
            %CorrectK1K2KBForPressure:
            deltaV    = nan(number_of_points,1); Kappa     = nan(number_of_points,1);
            lnK1fac   = nan(number_of_points,1); lnK2fac   = nan(number_of_points,1);
            lnKBfac   = nan(number_of_points,1);
            selected=(which_k1_k2==8);
            RR = (gas_constant.*temp_k);
            if any(selected)
                %***PressureEffectsOnK1inFreshWater:
                %               This is from Millero, 1983.
                deltaV(selected)  = -30.54 + 0.1849 .*temp_c(selected) - 0.0023366.*temp_c(selected).^2;
                Kappa(selected)   = (-6.22 + 0.1368 .*temp_c(selected) - 0.001233 .*temp_c(selected).^2)./1000;
                lnK1fac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
                %***PressureEffectsOnK2inFreshWater:
                %               This is from Millero, 1983.
                deltaV(selected)  = -29.81 + 0.115.*temp_c(selected) - 0.001816.*temp_c(selected).^2;
                Kappa(selected)   = (-5.74 + 0.093.*temp_c(selected) - 0.001896.*temp_c(selected).^2)./1000;
                lnK2fac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
                lnKBfac(selected) = 0 ;%; this doesn't matter since boron_concentration = 0 for this case
            end
            selected=(which_k1_k2==6 | which_k1_k2==7);
            if any(selected)
                %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
                %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
                %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
                %               but the fits are the same as those in
                %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
                %               who in turn quote Li, personal communication
                lnK1fac(selected) = (24.2 - 0.085.*temp_c(selected)).*pressure_bar(selected)./RR(selected);
                lnK2fac(selected) = (16.4 - 0.04 .*temp_c(selected)).*pressure_bar(selected)./RR(selected);
                %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
                %               and matches the GEOSECS results
                lnKBfac(selected) = (27.5 - 0.095.*temp_c(selected)).*pressure_bar(selected)./RR(selected);
            end
            selected=(which_k1_k2~=6 & which_k1_k2~=7 & which_k1_k2~=8);
            if any(selected)
                %***PressureEffectsOnK1:
                %               These are from Millero, 1995.
                %               They are the same as Millero, 1979 and Millero, 1992.
                %               They are from data of Culberson and Pytkowicz, 1968.
                deltaV(selected)  = -25.5 + 0.1271.*temp_c(selected);
                %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
                Kappa(selected)   = (-3.08 + 0.0877.*temp_c(selected))./1000;
                %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
	                lnK1fac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
                %               The fits given in Millero, 1983 are somewhat different.
                
                %***PressureEffectsOnK2:
                %               These are from Millero, 1995.
                %               They are the same as Millero, 1979 and Millero, 1992.
                %               They are from data of Culberson and Pytkowicz, 1968.
                deltaV(selected)  = -15.82 - 0.0219.*temp_c(selected);
                %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
                Kappa(selected)   = (1.13 - 0.1475.*temp_c(selected))./1000;
                %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
                lnK2fac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
                %               The fit given in Millero, 1983 is different.
                %               Not by a lot for deltaV, but by much for Kappa. %
                
                %***PressureEffectsOnKB:
                %               This is from Millero, 1979.
                %               It is from data of Culberson and Pytkowicz, 1968.
                deltaV(selected)  = -29.48 + 0.1622.*temp_c(selected) - 0.002608.*temp_c(selected).^2;
                %               Millero, 1983 has:
                %                 'deltaV = -28.56 + .1211.*temperature_in - .000321.*temperature_in.*temperature_in
                %               Millero, 1992 has:
                %                 'deltaV = -29.48 + .1622.*temperature_in + .295.*(Sali - 34.8)
                %               Millero, 1995 has:
                %                 'deltaV = -29.48 - .1622.*temperature_in - .002608.*temperature_in.*temperature_in
                %                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979
                Kappa(selected)   = -2.84./1000; % Millero, 1979
                %               Millero, 1992 and Millero, 1995 also have this.
                %                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
                %               Millero, 1983 has:
                %                 'Kappa = (-3 + .0427.*temperature_in)./1000
                lnKBfac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
            end
            k1_pressure_correction = exp(lnK1fac);
            k2_pressure_correction = exp(lnK2fac);
            kb_pressure_correction = exp(lnKBfac);
        end
        function kw_pressure_correction = calculate_pressure_correction_kw(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            RR = (gas_constant.*(temp_c+273.15));
            
            % CorrectKWForPressure:
            lnKWfac   = nan(number_of_points,1);
            deltaV = nan(number_of_points,1);
            Kappa = nan(number_of_points,1);
            selected=(which_k1_k2==8);
            if any(selected)
                % PressureEffectsOnKWinFreshWater:
                %               This is from Millero, 1983.
                deltaV(selected)  =  -25.6 + 0.2324.*temp_c(selected) - 0.0036246.*temp_c(selected).^2;
                Kappa(selected)   = (-7.33 + 0.1368.*temp_c(selected) - 0.001233 .*temp_c(selected).^2)./1000;
	                lnKWfac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
            
                %               NOTE the temperature dependence of KappaK1 and KappaKW
                %               for fresh water in Millero, 1983 are the same.
            end
            selected=(which_k1_k2~=8);
            if any(selected)
                % GEOSECS doesn't include OH term, so this won't matter.
                % Peng et al didn't include pressure, but here I assume that the KW correction
                %       is the same as for the other seawater cases.
                % PressureEffectsOnKW:
                %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
                deltaV(selected)  = -20.02 + 0.1119.*temp_c(selected) - 0.001409.*temp_c(selected).^2;
                %               Millero, 1992 and Millero, 1995 have:
                Kappa(selected)   = (-5.13 + 0.0794.*temp_c(selected))./1000; % Millero, 1983
                %               Millero, 1995 has this too, but Millero, 1992 is different.
                lnKWfac(selected) = (-deltaV(selected) + 0.5.*Kappa(selected).*pressure_bar(selected)).*pressure_bar(selected)./RR(selected);
                %               Millero, 1979 does not list values for these.
            end
            kw_pressure_correction = exp(lnKWfac);
        end
        function [kp1_pressure_correction,kp2_pressure_correction,kp3_pressure_correction] = calculate_pressure_correction_kp1_and_kp2_and_kp3(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            temp_k    = temp_c + 273.15;
                        
            % CorrectKP1KP2KP3KSiForPressure:
            % These corrections don't matter for the GEOSECS choice (which_k1_k2% = 6) and
            %       the freshwater choice (which_k1_k2% = 8). For the Peng choice I assume
            %       that they are the same as for the other choices (which_k1_k2% = 1 to 5).
            % The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
            %       same as Millero, 1983.
            % PressureEffectsOnKP1:
            deltaV = -14.51 + 0.1211.*temp_c - 0.000321.*temp_c.^2;
            Kappa  = (-2.67 + 0.0427.*temp_c)./1000;
            lnKP1fac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);
            % PressureEffectsOnKP2:
            deltaV = -23.12 + 0.1758.*temp_c - 0.002647.*temp_c.^2;
            Kappa  = (-5.15 + 0.09  .*temp_c)./1000;
            lnKP2fac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);
            % PressureEffectsOnKP3:
            deltaV = -26.57 + 0.202 .*temp_c - 0.003042.*temp_c.^2;
            Kappa  = (-4.08 + 0.0714.*temp_c)./1000;
            lnKP3fac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);

            kp1_pressure_correction = exp(lnKP1fac);
            kp2_pressure_correction = exp(lnKP2fac);
            kp3_pressure_correction = exp(lnKP3fac);
        end
        function ksi_pressure_correction = calculate_pressure_correction_ksi(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            temp_k    = temp_c + 273.15;
                        
            % PressureEffectsOnKSi:
            %  The only mention of this is Millero, 1995 where it is stated that the
            %    values have been estimated from the values of boric acid. HOWEVER,
            %    there is no listing of the values in the table.
            %    I used the values for boric acid from above.
            deltaV = -29.48 + 0.1622.*temp_c - 0.002608.*temp_c.^2;
            Kappa  = -2.84./1000;
            lnKSifac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);
            ksi_pressure_correction = exp(lnKSifac);
        end
        function knh4_pressure_correction = calculate_pressure_correction_knh4(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            temp_k    = temp_c + 273.15;
                        
            % PressureEffectsOnKNH4: added by J. Sharp
            deltaV = -26.43 + 0.0889.*temp_c - 0.000905.*temp_c.^2;
            Kappa  = (-5.03 + 0.0814.*temp_c)./1000;
            lnKNH4fac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);
            knh4_pressure_correction = exp(lnKNH4fac);
        end
        function kh2s_pressure_correction = calculate_pressure_correction_kh2s(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            temp_k    = temp_c + 273.15;

            % PressureEffectsOnKH2S: added by J. Sharp
            deltaV = -11.07 - 0.009.*temp_c - 0.000942.*temp_c.^2;
            Kappa  = (-2.89 + 0.054 .*temp_c)./1000;
            lnKH2Sfac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);
            kh2s_pressure_correction = exp(lnKH2Sfac);
        end

        function ks_pressure_correction = calculate_pressure_correction_ks(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            temp_k    = temp_c + 273.15;
                        
            % PressureEffectsOnKS:
            %       This is from Millero, 1995, which is the same as Millero, 1983.
            %       It is assumed that KS is on the free pH scale.
            deltaV = -18.03 + 0.0466.*temp_c + 0.000316.*temp_c.^2;
            Kappa = (-4.53 + 0.09.*temp_c)./1000;
            lnKSfac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);
            ks_pressure_correction = exp(lnKSfac);
        end
        function kf_pressure_correction = calculate_pressure_correction_kf(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar)
            temp_k    = temp_c + 273.15;
            
            % PressureEffectsOnKF:
            %       This is from Millero, 1995, which is the same as Millero, 1983.
            %       It is assumed that KF is on the free pH scale.
            deltaV = -9.78 - 0.009.*temp_c - 0.000942.*temp_c.^2;
            Kappa = (-3.91 + 0.054.*temp_c)./1000;
            lnKFfac = (-deltaV + 0.5.*Kappa.*pressure_bar).*pressure_bar./(gas_constant.*temp_k);

            kf_pressure_correction = exp(lnKFfac);
        end
        
        %% Combination
        function Ks = calculate_surface_all(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT)
            k0 = EquilibriumConstants.calculate_surface_k0(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
            [k1,k2] = EquilibriumConstants.calculate_surface_k1_and_k2(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
            kw = EquilibriumConstants.calculate_surface_kw(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
            kb = EquilibriumConstants.calculate_surface_kb(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
            [kp1,kp2,kp3,ksi] = EquilibriumConstants.calculate_surface_kp1_and_kp2_and_kp3_and_ksi(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
            [knh4,kh2s] = EquilibriumConstants.calculate_surface_knh4_and_kh2s(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
            
            ks = EquilibriumConstants.calculate_surface_ks(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2);
            kf = EquilibriumConstants.calculate_surface_kf(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2);
            
            Ks = containers.Map(["K0","K1","K2","KW","KB","KF","KS","KP1","KP2","KP3","KSi","KNH4","KH2S"], ...
                        {k0,k1,k2,kw,kb,kf,ks,kp1,kp2,kp3,ksi,knh4,kh2s});
        end
        function pressure_correction = calculate_pressure_correction_all(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar,co2_correction)
            k0_pressure_correction = EquilibriumConstants.calculate_pressure_correction_k0(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar,co2_correction);
            [k1_pressure_correction,k2_pressure_correction,kb_pressure_correction] = EquilibriumConstants.calculate_pressure_correction_k1_and_k2_and_kb(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);
            kw_pressure_correction = EquilibriumConstants.calculate_pressure_correction_kw(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);
            [kp1_pressure_correction,kp2_pressure_correction,kp3_pressure_correction] = EquilibriumConstants.calculate_pressure_correction_kp1_and_kp2_and_kp3(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);
            ksi_pressure_correction = EquilibriumConstants.calculate_pressure_correction_ksi(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);
            knh4_pressure_correction = EquilibriumConstants.calculate_pressure_correction_knh4(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);
            kh2s_pressure_correction = EquilibriumConstants.calculate_pressure_correction_kh2s(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);

            ks_pressure_correction = EquilibriumConstants.calculate_pressure_correction_ks(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);
            kf_pressure_correction = EquilibriumConstants.calculate_pressure_correction_kf(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar);

            pressure_correction = containers.Map(["K0","K1","K2","KW","KB","KF","KS","KP1","KP2","KP3","KSi","KNH4","KH2S"], ...
                        {k0_pressure_correction,k1_pressure_correction,k2_pressure_correction,kw_pressure_correction,kb_pressure_correction,kf_pressure_correction,ks_pressure_correction,kp1_pressure_correction,kp2_pressure_correction,kp3_pressure_correction,ksi_pressure_correction,knh4_pressure_correction,kh2s_pressure_correction});

        end
    
        function [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks)
            K0 = Ks("K0");
            K1 = Ks("K1");
            K2 = Ks("K2");
            KW = Ks("KW");
            KB = Ks("KB");
            KF = Ks("KF");
            KS = Ks("KS");
            KP1 = Ks("KP1");
            KP2 = Ks("KP2");
            KP3 = Ks("KP3");
            KSi = Ks("KSi");
            KNH4 = Ks("KNH4");
            KH2S = Ks("KH2S");
        end
        function Ks = pack_Ks(K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S)
                Ks = containers.Map(["K0","K1","K2","KW","KB","KF","KS","KP1","KP2","KP3","KSi","KNH4","KH2S"], ...
                        {K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S});
        end
    end
end