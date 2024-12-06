function [data,headers,nice_headers]=CO2SYS(parameter_1,parameter_2, ...
                                            parameter_1_type,parameter_2_type, ...
                                            salinity_in, ...
                                            temperature_in,temperature_out, ...
                                            pressure_in,pressure_out, ...
                                            silicate,phosphate,ammonia,sulphide, ...
                                            pH_scale_in, ...
                                            which_k1_k2,which_kso4,which_kf, which_boron, ...
                                            varargin)
    
    % Input conditioning

    % Determine lengths of input vectors
    input_lengths = [length(parameter_1) length(parameter_2) length(parameter_1_type)...
                length(parameter_2_type) length(salinity_in) length(temperature_in)...
                length(temperature_out) length(pressure_in) length(pressure_out)...
                length(silicate) length(phosphate) length(ammonia) length(sulphide)...
                length(pH_scale_in) length(which_k1_k2) length(which_kso4)...
	            length(which_kf) length(which_boron)];
    
    if length(unique(input_lengths))>2
	    disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end
    
    % Find the longest column vector:
    number_of_points = max(input_lengths);
        
    % set default for optional input argument
    co2_pressure_correction = zeros(numel(number_of_points),1);
    % parse optional input argument
    for index = 1:2:length(varargin)-1
        if strcmpi(varargin{index},'co2_press')
            co2_pressure_correction = varargin{index+1};
        end
    end
    
    % Populate column vectors
    parameter_1(1:number_of_points,1)      = parameter_1(:);
    parameter_2(1:number_of_points,1)      = parameter_2(:);
    parameter_1_type(1:number_of_points,1) = parameter_1_type(:);
    parameter_2_type(1:number_of_points,1) = parameter_2_type(:);
    salinity_in(1:number_of_points,1)      = salinity_in(:);
    temperature_in(1:number_of_points,1)   = temperature_in(:);
    temperature_out(1:number_of_points,1)  = temperature_out(:);
    pressure_in(1:number_of_points,1)      = pressure_in(:);
    pressure_out(1:number_of_points,1)     = pressure_out(:) ;
    silicate(1:number_of_points,1)         = silicate(:);
    phosphate(1:number_of_points,1)        = phosphate(:);
    ammonia(1:number_of_points,1)          = ammonia(:);
    sulphide(1:number_of_points,1)         = sulphide(:);
    pH_scale_in(1:number_of_points,1)      = pH_scale_in(:);
    which_k1_k2(1:number_of_points,1)      = which_k1_k2(:);
    which_kso4(1:number_of_points,1)       = which_kso4(:);
    which_kf(1:number_of_points,1)         = which_kf(:);
    which_boron(1:number_of_points,1)      = which_boron(:);
    co2_pressure_correction(1:number_of_points,1) = co2_pressure_correction(:);
    
    % Generate empty vectors for...
    alkalinity = NaN(number_of_points,1);
    dic        = NaN(number_of_points,1);
    pH         = NaN(number_of_points,1);
    pco2       = NaN(number_of_points,1);
    fco2       = NaN(number_of_points,1);
    hco3       = NaN(number_of_points,1);
    co3        = NaN(number_of_points,1);
    co2        = NaN(number_of_points,1);

    which_ks = WhichKs(which_k1_k2,which_kso4,which_kf);

    salinity = salinity_in;
    
    % Assign values to empty vectors.
    selected = (parameter_1_type==1 & parameter_1~=-999);   
    alkalinity(selected) = parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg

    selected = (parameter_1_type==2 & parameter_1~=-999);   
    dic(selected) = parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg

    selected = (parameter_1_type==3 & parameter_1~=-999);   
    pH(selected) = parameter_1(selected);

    selected = (parameter_1_type==4 & parameter_1~=-999);   
    pco2(selected) = parameter_1(selected)/1e6; % Convert from microatm. to atm.

    selected = (parameter_1_type==5 & parameter_1~=-999);   
    fco2(selected) = parameter_1(selected)/1e6; % Convert from microatm. to atm.

    selected = (parameter_1_type==6 & parameter_1~=-999); 
    hco3(selected) = parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_1_type==7 & parameter_1~=-999);  
    co3(selected) = parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_1_type==8 & parameter_1~=-999);  
    co2(selected) = parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_2_type==1 & parameter_2~=-999);   
    alkalinity(selected) = parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_2_type==2 & parameter_2~=-999);   
    dic(selected) = parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_2_type==3 & parameter_2~=-999);   
    pH(selected) = parameter_2(selected);
    
    selected = (parameter_2_type==4 & parameter_2~=-999);   
    pco2(selected) = parameter_2(selected)/1e6; % Convert from microatm. to atm.
    
    selected = (parameter_2_type==5 & parameter_2~=-999);   
    fco2(selected) = parameter_2(selected)/1e6; % Convert from microatm. to atm.
    
    selected = (parameter_2_type==6 & parameter_2~=-999); 
    hco3(selected) = parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_2_type==7 & parameter_2~=-999);  
    co3(selected) = parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected = (parameter_2_type==8 & parameter_2~=-999);  
    co2(selected) = parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    % Generate the columns holding Si, Phos, Amm, H2S and salinity.
    % Pure Water case:
    selected = (which_k1_k2==8);
    salinity(selected) = 0;

    composition = Composition(salinity)...
                    .set_silicate_concentration(silicate/1e6)...
                    .set_phosphate_concentration(phosphate/1e6)...
                    .set_ammonia_concentration(ammonia/1e6)...
                    .set_sulphide_concentration(sulphide/1e6)...
                    .estimate_all_from_salinity(which_boron)...
                    .remove_freshwater_species(which_ks)...
                    .adjust_geosecs_species(which_ks)...
                    .calculate_peng_correction(which_ks);

    % Calculate the constants for all samples at input conditions
    Ks_in = Ks(temperature_in,...
               pressure_in/10,...
               composition,...
               pH_scale_in,...
               which_ks,...
               co2_pressure_correction)...
             .calculate_all();

    K0 = Ks_in.k0;
    
    % Make sure fCO2 is available for each sample that has pCO2 or CO2.
    temp_k = temperature_in+273.15;
    fugacity_factor = calculate_fugacity_factor(co2_pressure_correction,number_of_points,which_k1_k2,temp_k);
    
    selected = (~isnan(pco2) & (parameter_1_type==4 | parameter_2_type==4));
    fco2(selected) = pco2(selected).*fugacity_factor(selected);

    selected = (~isnan(co2) & (parameter_1_type==8 | parameter_2_type==8)); 
    fco2(selected) = co2(selected)./K0(selected);
    
    % Generate vectors for results, and copy the raw input values into them
    alkalinity_in = alkalinity;
    dic_in        = dic;
    pH_in         = pH;
    pco2_in       = pco2;
    fco2_in       = fco2;
    hco3_in       = hco3;
    co3_in        = co3;
    co2_in        = co2;
    
    % Generate vector describing the combination of input parameters
    % So, the valid ones are:
    % 12,13,15,16,17,18,23,25,26,27,28,35,36,37,38,56,57,67,68,78
    combination = 10*min(parameter_1_type,parameter_2_type) + max(parameter_1_type,parameter_2_type);
    
    % Calculate missing values for AT,CT,PH,FC,HCO3,CO3,CO2:
    % pCO2 will be calculated later on, routines work with fCO2.
    selected = (combination==12); % input TA, TC
    if any(selected)
    selected = (~isnan(alkalinity_in) & ~isnan(dic_in) & selected);
        pH_in(selected) = calculate_pH_from_alkalinity_dic(alkalinity_in(selected)-composition.peng_correction(selected),dic_in(selected),Ks_in.select(selected));

        selected = (~isnan(pH_in) & selected);
        if any(selected)
           fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
           hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
           co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        end
    end

    selected = (combination==13); % input TA, pH
    if any(selected)
        selected = (~isnan(alkalinity_in) & ~isnan(pH_in) & selected);
        dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
    end

    selected = (combination==14 | combination==15 | combination==18); % input TA, (pCO2 or fCO2 or CO2)
    if any(selected)
        selected=(~isnan(alkalinity_in) & ~isnan(fco2_in) & selected);
        pH_in(selected) = calculate_pH_from_alkalinity_fco2(alkalinity_in(selected)-composition.peng_correction(selected),fco2_in(selected),Ks_in.select(selected));
        selected = (~isnan(pH_in) & selected);
        if any(selected)
            dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
            hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
            co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        end
    end

    selected = (combination==16); % input TA, HCO3
    if any(selected)
    selected = (~isnan(alkalinity_in) & ~isnan(hco3_in) & selected);
        pH_in(selected) = calculate_pH_from_alkalinity_hco3(alkalinity_in(selected)-composition.peng_correction(selected),hco3_in(selected),Ks_in.select(selected));  % added Peng correction // MPH
        selected=(~isnan(pH_in) & selected);
        if any(selected)
           dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
           fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
           co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
        end
    end

    selected = (combination==17); % input TA, CO3
    if any(selected)
    selected = (~isnan(alkalinity_in) & ~isnan(co3_in) & selected);
        pH_in(selected) = calculate_pH_from_alkalinity_co3(alkalinity_in(selected)-composition.peng_correction(selected),co3_in(selected),Ks_in.select(selected));  % added Peng correction // MPH
        selected=(~isnan(pH_in) & selected);
        if any(selected)
           dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
           fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
           hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
        end
    end

    selected = (combination==23); % input TC, pH
    if any(selected)
        selected = (~isnan(dic_in) & ~isnan(pH_in) & selected);
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
    end
    
    selected = (combination==24 | combination==25 | combination==28);  % input TC, (pCO2 or fCO2 or CO2)
    if any(selected)
        selected = (~isnan(dic_in) & ~isnan(fco2_in) & selected);
        pH_in(selected) = calculate_pH_from_dic_fco2(dic_in(selected),fco2_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
    end

    selected = (combination==26); % input TC, HCO3
    if any(selected)
    selected = (~isnan(dic_in) & ~isnan(hco3_in) & selected);
        pH_in(selected) = calculate_pH_from_dic_hco3(dic_in(selected), hco3_in(selected), Ks_in.select(selected)); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected), pH_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
    end

    selected = (combination==27); % input TC, CO3
    if any(selected)
    selected = (~isnan(dic_in) & ~isnan(co3_in) & selected);
        % [pH_in(selected),fco2_in(selected)] = calculate_pH_fco2_from_dic_co3(dic_in(selected),co3_in(selected), Ks_in.select(selected));
        pH_in(selected) = calculate_pH_from_dic_co3(dic_in(selected), co3_in(selected), Ks_in.select(selected)); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected), pH_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
    end

    selected = (combination==34 | combination==35 | combination==38); % input pH, (pCO2 or fCO2 or CO2)
    if any(selected)
        selected = (~isnan(pH_in) & ~isnan(fco2_in) & selected);
        dic_in(selected) = calculate_dic_from_pH_fco2(pH_in(selected),fco2_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
    end

    selected = (combination==36); % input pH, HCO3
    if any(selected)
    selected = (~isnan(pH_in) & ~isnan(hco3_in) & selected);
        alkalinity_in(selected) = calculate_alkalinity_from_pH_hco3(pH_in(selected),hco3_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
    end

    selected = (combination==37); % input pH, CO3
    if any(selected)
    selected = (~isnan(pH_in) & ~isnan(co3_in) & selected);
        alkalinity_in(selected) = calculate_alkalinity_from_pH_co3(pH_in(selected),co3_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
    end

    selected = (combination==46 | combination==56 | combination==68); % input (pCO2 or fCO2 or CO2), HCO3
    if any(selected)
    selected = (~isnan(fco2_in) & ~isnan(hco3_in) & selected);
        pH_in(selected) = calculate_pH_from_fco2_hco3(fco2_in(selected),hco3_in(selected), Ks_in.select(selected));
        dic_in(selected) = calculate_dic_from_pH_fco2(pH_in(selected),fco2_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        co3_in(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
    end

    selected = (combination==47 | combination==57 | combination==78); % input (pCO2 or fCO2 or CO2), CO3
    if any(selected)
    selected = (~isnan(fco2_in) & ~isnan(co3_in) & selected);
        pH_in(selected) = calculate_pH_from_fco2_co3(fco2_in(selected),co3_in(selected), Ks_in.select(selected));
        dic_in(selected) = calculate_dic_from_pH_fco2 (pH_in(selected),fco2_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        hco3_in(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_in(selected), Ks_in.select(selected));
    end

    selected = (combination==67); % input HCO3, CO3
    if any(selected)
    selected = (~isnan(hco3_in) & ~isnan(co3_in) & selected);
        pH_in(selected) = calculate_pH_from_co3_hco3(co3_in(selected),hco3_in(selected), Ks_in.select(selected));
        alkalinity_in(selected) = calculate_alkalinity_from_pH_co3(pH_in(selected),co3_in(selected),Ks_in.select(selected)) + composition.peng_correction(selected);
        dic_in(selected) = calculate_dic_from_alkalinity_pH(alkalinity_in(selected)-composition.peng_correction(selected),pH_in(selected),Ks_in.select(selected));
        fco2_in(selected) = calculate_fco2_from_dic_pH(dic_in(selected),pH_in(selected),Ks_in.select(selected));
        %CO2ic(selected)                = CalculateCO2fromTCpH(TCc(selected),PHic(selected));
    end
    
    % By now, an fCO2 value is available for each sample.
    % Generate the associated pCO2 value:
    selected = (isnan(pco2_in) & (parameter_1_type~=4 | parameter_2_type~=4)); 
    pco2_in(selected)  = fco2_in(selected)./fugacity_factor(selected);

    % Generate the associated CO2 value:
    selected = (isnan(co2_in) & (parameter_1_type~=8 | parameter_2_type~=8)); 
    co2_in(selected) = fco2_in(selected).*K0(selected);
    
    % Calculate Other Params At Input Conditions
    % Generate vector of NaN
    nan_vector = NaN(number_of_points,1);
    % Copy that vector into the alkalinity components
    [boron_alkalinity_in,oh_alkalinity_in,phosphate_alkalinity_in,...
        silicate_alkalinity_in,ammonia_alkalinity_in,sulphide_alkalinity_in,...
        h_free_alkalinity_in,sulphate_alkalinity_in,fluorine_alkalinity_in,...
        revelle_alkalinity_in,saturation_state_calcite_in,saturation_state_aragonite_in,...
        co2_dry_alkalinity_in] = deal(nan_vector);

    selected = (~isnan(pH_in)); % if PHic = NaN, pH calculation was not performed or did not converge
    [boron_alkalinity_in(selected),oh_alkalinity_in(selected),phosphate_alkalinity_in(selected),...
        silicate_alkalinity_in(selected),ammonia_alkalinity_in(selected),sulphide_alkalinity_in(selected),...
        h_free_alkalinity_in(selected),sulphate_alkalinity_in(selected),fluorine_alkalinity_in(selected),...
    ] = calculate_alkalinity_parts(pH_in(selected),Ks_in.select(selected));
    
    phosphate_alkalinity_in(selected) = phosphate_alkalinity_in(selected)+composition.peng_correction(selected);
    revelle_alkalinity_in(selected) = calculate_revelle_factor(alkalinity_in(selected)-composition.peng_correction(selected), dic_in(selected),Ks_in.select(selected));
    [saturation_state_calcite_in(selected),saturation_state_aragonite_in(selected)] = calculate_carbonate_solubility(salinity(selected), temperature_in(selected), dic_in(selected), pH_in(selected), Ks_in.select(selected),composition.calcium(selected),which_k1_k2(selected),pressure_in(selected)/10);
    vapour_pressure_factor = calculate_vapour_pressure_factor(salinity,temp_k);
    co2_dry_alkalinity_in(~isnan(pco2_in),1) = pco2_in(~isnan(pco2_in),1)./vapour_pressure_factor(~isnan(pco2_in),1); % ' this assumes pTot = 1 atm
    
    
    substrate_inhibitor_ratio_in = hco3_in./(h_free_alkalinity_in.*1e6);
    
    % % Just for reference, convert pH at input conditions to the other scales
    pHicT = NaN(number_of_points,1);
    pHicS = NaN(number_of_points,1);
    pHicF = NaN(number_of_points,1);
    pHicN = NaN(number_of_points,1);
    relevant_ks = Ks_in.select(selected);
    [pHicT(selected),pHicS(selected),pHicF(selected),pHicN(selected)]=relevant_ks.controls.pH_scale_conversion(2).find_pH_on_all_scales(pH_in(selected),relevant_ks.controls);
    
    % Merge the Ks at input into an array. Ks at output will be glued to this later.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KC,KNH4,KH2S] = Ks_in.unpack_all();
    k_in_vector = [K0,K1,K2,-log10(K1),-log10(K2),KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S];

    clear K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KC KNH4 KH2S
    
    % Calculate the constants for all samples at output conditions
    Ks_out = Ks(temperature_out,...
                pressure_out/10,...
                composition,...
                pH_scale_in,...
                which_ks,...
                co2_pressure_correction)...
              .calculate_all();


    % For output conditions, using conservative TA and TC, calculate pH, fCO2
    % and pCO2, HCO3, CO3, and CO2
    temp_k = temperature_out+273.15;
    selected=(~isnan(alkalinity_in) & ~isnan(dic_in)); % i.e., do for all samples that have TA and TC values
    pH_out = NaN(number_of_points,1);
    [co3_out,hco3_out,fco2_out] = deal(pH_out);
    pH_out(selected) = calculate_pH_from_alkalinity_dic(alkalinity_in(selected)-composition.peng_correction(selected), dic_in(selected),Ks_out.select(selected)); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
    fco2_out(selected) = calculate_fco2_from_dic_pH(dic_in(selected), pH_out(selected), Ks_out.select(selected));
    hco3_out(selected) = calculate_hco3_from_dic_pH(dic_in(selected),pH_out(selected),Ks_out.select(selected));
    co3_out(selected) = calculate_co3_from_dic_pH(dic_in(selected),pH_out(selected),Ks_out.select(selected));
    
    % Generate the associated pCO2 value:
    fugacity_factor = calculate_fugacity_factor(co2_pressure_correction,number_of_points,which_k1_k2,temp_k);
    pco2_out  = fco2_out./fugacity_factor;
    % Generate the associated CO2 value:

    K0_out = Ks_out.k0;
    co2_out = fco2_out.*K0_out;
    
    % Calculate Other Params At Output Conditions
    [boron_alkalinity_out,oh_alkalinity_out,phosphate_alkalinity_out,...
        silicate_alkalinity_out,ammonia_alkalinity_out,sulphide_alkalinity_out,...
        h_free_alkalinity_out,sulphate_alkalinity_out,fluorine_alkalinity_out,...
        revelle_alkalinity_out,saturation_state_calcite_out,saturation_state_aragonite_out,...
        co2_dry_alkalinity_out] = deal(nan_vector);

    [boron_alkalinity_out(selected),oh_alkalinity_out(selected),phosphate_alkalinity_out(selected),silicate_alkalinity_out(selected),ammonia_alkalinity_out(selected),...
        sulphide_alkalinity_out(selected), h_free_alkalinity_out(selected),sulphate_alkalinity_out(selected),fluorine_alkalinity_out(selected)] = calculate_alkalinity_parts(pH_out(selected),Ks_out.select(selected));
    
    phosphate_alkalinity_out(selected)                 = phosphate_alkalinity_out(selected)+composition.peng_correction(selected);
    revelle_alkalinity_out(selected)              = calculate_revelle_factor(alkalinity_in(selected)-composition.peng_correction(selected), dic_in(selected),Ks_out.select(selected));
    [saturation_state_calcite_out(selected),saturation_state_aragonite_out(selected)] = calculate_carbonate_solubility(salinity(selected), temperature_out(selected), dic_in(selected), pH_out(selected), Ks_out.select(selected), composition.calcium(selected),which_k1_k2(selected),pressure_out(selected)/10);
    vapour_pressure_factor = calculate_vapour_pressure_factor(salinity,temp_k);
    co2_dry_alkalinity_out(~isnan(pco2_out),1)    = pco2_out(~isnan(pco2_out))./vapour_pressure_factor(~isnan(pco2_out)); % ' this assumes pTot = 1 atm
    substrate_inhibitor_ratio_out = hco3_out./(h_free_alkalinity_out.*1e6);
    
    % Just for reference, convert pH at output conditions to the other scales
    pH_out_total = NaN(number_of_points,1);
    pH_out_seawater = NaN(number_of_points,1);
    pH_out_free = NaN(number_of_points,1);
    pH_out_NBS = NaN(number_of_points,1);
        
    relevant_ks = Ks_out.select(selected);
    [pH_out_total(selected),pH_out_seawater(selected),pH_out_free(selected),pH_out_NBS(selected)] = relevant_ks.controls.pH_scale_conversion(2).find_pH_on_all_scales(pH_out(selected),relevant_ks.controls);

    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KC,KNH4,KH2S] = Ks_out.unpack_all();
    k_out_vector = [K0,K1,K2,-log10(K1),-log10(K2),KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S];
    concentration_vector =[composition.boron,composition.fluorine,composition.sulphate,composition.phosphate,composition.silicate,composition.ammonia,composition.sulphide];
    
    % Saving data in array, 99 columns, as many rows as samples input
    data=[alkalinity_in*1e6         dic_in*1e6        pH_in           pco2_in*1e6        fco2_in*1e6...
          hco3_in*1e6      co3_in*1e6      co2_in*1e6      boron_alkalinity_in*1e6     oh_alkalinity_in*1e6...
          phosphate_alkalinity_in*1e6     silicate_alkalinity_in*1e6   ammonia_alkalinity_in*1e6  sulphide_alkalinity_in*1e6    h_free_alkalinity_in*1e6... %%% Multiplied Hfreeinp *1e6, svh20100827
          revelle_alkalinity_in      saturation_state_calcite_in     saturation_state_aragonite_in     co2_dry_alkalinity_in*1e6  substrate_inhibitor_ratio_in...
          pH_out            pco2_out*1e6       fco2_out*1e6       hco3_out*1e6      co3_out*1e6...
          co2_out*1e6       boron_alkalinity_out*1e6    oh_alkalinity_out*1e6      phosphate_alkalinity_out*1e6     silicate_alkalinity_out*1e6...
          ammonia_alkalinity_out*1e6   sulphide_alkalinity_out*1e6   h_free_alkalinity_out*1e6   revelle_alkalinity_out      saturation_state_calcite_out... %%% Multiplied Hfreeout *1e6, svh20100827
          saturation_state_aragonite_out      co2_dry_alkalinity_out*1e6 substrate_inhibitor_ratio_out         pHicT           pHicS...
          pHicF           pHicN          pH_out_total          pH_out_seawater           pH_out_free...
          pH_out_NBS           temperature_in         temperature_out        pressure_in          pressure_out...
          parameter_1_type        parameter_2_type       which_k1_k2  which_kso4    which_kf...
          which_boron           pH_scale_in      salinity_in            phosphate             silicate...
          ammonia             sulphide            k_in_vector          k_out_vector           concentration_vector*1e6];
    data(isnan(data))=-999;
    
    headers={'TAlk';'TCO2';'pHin';'pCO2in';'fCO2in';'HCO3in';'CO3in';...
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
        'KP3output';'KSioutput';'KNH4output';'KH2Soutput';'boron_concentration';'fluorine_concentration';'sulphate_concentration';...
        'phosphate_concentration';'silicate_concentration';'ammonia_concentration';'sulphide_concentration'};
    
    nice_headers={...
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
        '93 - boron_concentration               (umol/kgSW) ';
        '94 - fluorine_concentration               (umol/kgSW) ';
        '95 - sulphate_concentration               (umol/kgSW) ';
        '96 - phosphate_concentration               (umol/kgSW) ';
        '97 - silicate_concentration              (umol/kgSW) ';
        '98 - ammonia_concentration             (umol/kgSW) ';
        '99 - sulphide_concentration             (umol/kgSW) '};	
end 

%% Calculate pH
function pH_out = calculate_pH_from_alkalinity_dic(alkalinity,dic,Ks)
    [K1,K2,KB] = Ks.unpack_some(["K1","K2","KB"]);
    
    % Find initital pH guess using method of Munhoven (2013)
    initial_pH_guess = calculate_pH_from_alkalinity_dic_munhoven(alkalinity, dic, Ks);
    pH = initial_pH_guess;
    pH_tolerance = 1e-4;  % tolerance for iterations end
    
    counter = 0;
    delta_pH(1:numel(alkalinity),1) = pH_tolerance+1;
    above_tolerance = (abs(delta_pH) > pH_tolerance);

    
    while any(above_tolerance)
        H = 10.^(-pH);
        carbonate_denominator = (H.^2 + K1.*H + K1.*K2);
        carbonate_alkalinity = dic.*K1.*(H + 2.*K2)./carbonate_denominator;        
        
        [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

        residual  = alkalinity - carbonate_alkalinity - boron_alkalinity - hydroxide_alkalinity - phosphate_alkalinity - silicate_alkalinity  - ammonia_alkalinity - sulphide_alkalinity + hydrogen_free + sulphate_acidity + fluorine_acidity;
        
        % find Slope dTA/dpH;
        % (this is not exact, but keeps all important terms);
        slope = log(10).*(dic.*K1.*H.*(H.*H + K1.*K2 + 4.*H.*K2)./carbonate_denominator./carbonate_denominator + boron_alkalinity.*H./(KB + H) + hydroxide_alkalinity + H);
        delta_pH = residual./slope; %' this is Newton's method
        
        % ' to keep the jump from being too big:
        while any(abs(delta_pH) > 1)
            high_delta = abs(delta_pH)>1; 
            delta_pH(high_delta) = delta_pH(high_delta)/2;
        end

        pH(above_tolerance) = pH(above_tolerance) + delta_pH(above_tolerance);
        above_tolerance     = abs(delta_pH) > pH_tolerance;
        counter=counter+1;
     
        if counter>10000
            failed = find(abs(delta_pH) > pH_tolerance);
            pH(failed) = NaN;
            show_failed(failed);
        end
    end
    pH_out = pH;
end

function pH_out = calculate_pH_from_alkalinity_fco2(alkalinity, fco2,Ks)
    [K0,K1,K2,KB] = Ks.unpack_some(["K0","K1","K2","KB"]);

    % Find initital pH guess using method of Munhoven (2013)
    co2 = fco2.*K0; % Convert fCO2 to CO2
    
    pH_initial_guess = calculate_pH_from_alkalinity_co2_munhoven(alkalinity, co2, Ks);
    pH = pH_initial_guess;
    pH_tolerance = 1e-4;

    delta_pH = pH_tolerance + 1;
    counter = 0;
    

    above_tolerance = (abs(delta_pH) > pH_tolerance);
    while any(above_tolerance)
        H = 10.^(-pH);
        hco3 = K0.*K1.*fco2./H;
        co3 = K0.*K1.*K2.*fco2./(H.*H);
        carbonate_alkalinity = hco3 + 2.*co3;

        [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);
        
        residual  = alkalinity - carbonate_alkalinity - boron_alkalinity - hydroxide_alkalinity - phosphate_alkalinity - silicate_alkalinity - ammonia_alkalinity - sulphide_alkalinity + hydrogen_free + sulphate_acidity + fluorine_acidity;
        slope = log(10).*(hco3 + 4.*co3 + boron_alkalinity.*H./(KB + H) + hydroxide_alkalinity + H);
        delta_pH   = residual./slope;

        while any(abs(delta_pH) > 1)
            high_delta = abs(delta_pH)>1; 
            delta_pH(high_delta) = delta_pH(high_delta)/2;
        end

        pH(above_tolerance) = pH(above_tolerance) + delta_pH(above_tolerance);
        above_tolerance = abs(delta_pH) > pH_tolerance;
        counter = counter+1;
     
        if counter>10000
            failed = find(abs(delta_pH) > pH_tolerance);
            pH(failed) = NaN;
            show_failed(failed);
        end
    end
    pH_out = pH;
end

function pH_out = calculate_pH_from_dic_fco2(dic,fco2,Ks)
    [K0,K1,K2] = Ks.unpack_some(["K0","K1","K2"]);
    RR = K0.*fco2./dic;
    %       if RR >= 1
    %          varargout{1}= missingn;
    %          disp('nein!');return;
    %       end
    % check after sub to see if pH = missingn.
    Discr = (K1.*RR).*(K1.*RR) + 4.*(1 - RR).*(K1.*K2.*RR);
    H     = 0.5.*(K1.*RR + sqrt(Discr))./(1 - RR);
    %       if (H <= 0)
    %           pHctemp = missingn;
    %       else
    pH_out = log(H)./log(0.1);
    %       end
end

function pH_out = calculate_pH_from_alkalinity_hco3(alkalinity, hco3,Ks)
    [K2,KB] = Ks.unpack_some(["K2","KB"]);

    pH_initial_guess = calculate_pH_from_alkalinity_hco3_munhoven(alkalinity,hco3,Ks);
    pH = pH_initial_guess;
    pH_tolerance = 1e-4; % tolerance
    delta_pH = pH_tolerance+1;

    counter = 0;


    above_tolerance = (abs(delta_pH) > pH_tolerance);
    while any(above_tolerance)
        H = 10.^(-pH);
        carbonate_alkalinity = hco3.*(H+2.*K2)./H;

        [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

        residual  = alkalinity - carbonate_alkalinity - boron_alkalinity - hydroxide_alkalinity - phosphate_alkalinity - silicate_alkalinity - ammonia_alkalinity - sulphide_alkalinity + hydrogen_free + sulphate_acidity + fluorine_acidity;

        slope = log(10) .* (2 .* hco3 .* K2 ./ H + boron_alkalinity .* H ./ (KB + H) + hydroxide_alkalinity + H);
        
        delta_pH   = residual./slope; %' this is Newton's method
        
        while any(abs(delta_pH) > 1)
            high_delta = abs(delta_pH)>1;
            delta_pH(high_delta) = delta_pH(high_delta)./2;
        end
        pH(above_tolerance) = pH(above_tolerance) + delta_pH(above_tolerance);
        above_tolerance = abs(delta_pH) > pH_tolerance;

        counter = counter+1;
     
        if counter>10000
            failed = find(abs(delta_pH) > pH_tolerance);
            pH(failed) = NaN;
            show_failed(failed);
        end
    end
    pH_out = pH;
end

function pH_out = calculate_pH_from_dic_hco3(dic, hco3, Ks)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);
    
    RR = dic./hco3;
    discriminant = ((1-RR).*(1-RR) - 4.*(1./(K1)).*(K2));
    H = 0.5.*((-(1-RR)) - sqrt(discriminant))./(1./(K1)); % Subtraction
    
    pH_out = log(H)./log(0.1);
end

function pH_out = calculate_pH_from_alkalinity_co3(alkalinity,co3,Ks)
    [K2,KB] = Ks.unpack_some(["K2","KB"]);

    pH_initial_guess = calculate_pH_from_alkalinity_co3_munhoven(alkalinity,co3,Ks);
    pH = pH_initial_guess;
    pH_tolerance = 1e-4; % tolerance
    delta_pH = pH_tolerance + 1.0;

    counter = 0;

    
    above_tolerance = (abs(delta_pH) > pH_tolerance);
    while any(above_tolerance)
        H = 10.^(-pH);
        carbonate_alkalinity = co3.*(H+2.*K2)./K2;
        
        [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

        residual  = alkalinity - carbonate_alkalinity - boron_alkalinity - hydroxide_alkalinity - phosphate_alkalinity - silicate_alkalinity - ammonia_alkalinity - sulphide_alkalinity + hydrogen_free + sulphate_acidity + fluorine_acidity;


        slope = log(10) .* (-co3 .* H ./ K2 + boron_alkalinity .* H ./ (KB + H) + hydroxide_alkalinity + H);
        delta_pH   = residual./slope;

        while any(abs(delta_pH) > 1)
            high_delta = abs(delta_pH)>1; 
            delta_pH(high_delta) = delta_pH(high_delta)/2;
        end
        
        pH(above_tolerance) = pH(above_tolerance) + delta_pH(above_tolerance);
        above_tolerance     = abs(delta_pH) > pH_tolerance;
        counter=counter+1;
     
        if counter>10000 
            failed = find(abs(delta_pH) > pH_tolerance);
            pH(failed) = NaN;
            show_failed(failed);
        end
    end
    pH_out = pH;
end

function pH_out = calculate_pH_from_dic_co3(dic, co3, Ks)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);

    RR = dic./co3;

    discriminant = ((1./K2).*(1./K2) - 4.*(1./(K1.*K2)).*(1-RR));
    H = 0.5.*((-1./K2) + sqrt(discriminant))./(1./(K1.*K2)); % Addition
    pH_out = log(H)./log(0.1);
end

function pH_out = calculate_pH_from_fco2_co3(fco2, co3, Ks)
    [K0,K1,K2] = Ks.unpack_some(["K0","K1","K2"]);
    H = sqrt((fco2.*K0.*K1.*K2)./co3);    % removed incorrect (selected) index from CO3i // MPH
    pH_out = -log10(H);
end

function pH_out = calculate_pH_from_fco2_hco3(fco2, hco3, Ks)
    [K0,K1] = Ks.unpack_some(["K0","K1"]);
    H = (fco2.*K0.*K1)./hco3;  % removed incorrect (selected) index from HCO3i // MPH
    pH_out = -log10(H);
end

function pH_out = calculate_pH_from_co3_hco3(co3, hco3, Ks)
    K2 = Ks.k2;
    H = hco3.*K2./co3;
    pH_out = -log10(H);
end


%% Calculate alkalinity
function alkalinity = calculate_alkalinity_from_dic_pH(dic,pH,Ks)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);

    H = 10.^(-pH);
    carbonate_alkalinity = dic.*K1.*(H + 2.*K2)./(H.*H + K1.*H + K1.*K2);

    [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

    
    alkalinity = carbonate_alkalinity + boron_alkalinity + hydroxide_alkalinity + phosphate_alkalinity + silicate_alkalinity + ammonia_alkalinity + sulphide_alkalinity - hydrogen_free - sulphate_acidity - fluorine_acidity;
end

function alkalinity = calculate_alkalinity_from_pH_hco3(pH,hco3,Ks)
    K2 = Ks.k2;

    H = 10.^(-pH);
    carbonate_alkalinity = hco3.*(2.*K2./H + 1);    
    
    [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

    alkalinity = carbonate_alkalinity + boron_alkalinity + hydroxide_alkalinity + phosphate_alkalinity + silicate_alkalinity + ammonia_alkalinity + sulphide_alkalinity - hydrogen_free - sulphate_acidity - fluorine_acidity;
end

function alkalinity = calculate_alkalinity_from_pH_co3(pH,co3,Ks)
    K2 = Ks.k2;

    H = 10.^(-pH);
    carbonate_alkalinity = co3.*(H./K2 + 2);    
    
    [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

    
    alkalinity = carbonate_alkalinity + boron_alkalinity + hydroxide_alkalinity + phosphate_alkalinity + silicate_alkalinity + ammonia_alkalinity + sulphide_alkalinity - hydrogen_free - sulphate_acidity - fluorine_acidity;
end


%% Calculate dic
function dic = calculate_dic_from_alkalinity_pH(alkalinity,pH,Ks)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);
    [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks);

    H = 10.^(-pH);

    carbonate_alkalinity = alkalinity - boron_alkalinity - hydroxide_alkalinity - phosphate_alkalinity - silicate_alkalinity - ammonia_alkalinity - sulphide_alkalinity + hydrogen_free + sulphate_acidity + fluorine_acidity;
    dic = carbonate_alkalinity.*(H.*H + K1.*H + K1.*K2)./(K1.*(H + 2.*K2));
end

function dic = calculate_dic_from_pH_fco2(pH,fco2,Ks)
    [K0,K1,K2] = Ks.unpack_some(["K0","K1","K2"]);

    H = 10.^(-pH);
    dic = K0.*fco2.*(H.*H + K1.*H + K1.*K2)./(H.*H);
end


%% Calculate carbon species
function fco2 = calculate_fco2_from_dic_pH(dic,pH,Ks)
    [K0,K1,K2] = Ks.unpack_some(["K0","K1","K2"]);

    H = 10.^(-pH);
    fco2 = dic.*H.*H./(H.*H + K1.*H + K1.*K2)./K0;
end

function co3 = calculate_co3_from_dic_pH(dic, pH, Ks)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);
    H = 10.^(-pH);
    co3 = dic.*K1.*K2./(K1.*H + H.*H + K1.*K2);
end

function hco3 = calculate_hco3_from_dic_pH(dic, pH, Ks)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);
    H = 10.^(-pH);
    hco3 = dic.*K1.*H./(K1.*H + H.*H + K1.*K2);
end

%% Munhovens
function pH_out = calculate_pH_from_alkalinity_dic_munhoven(alkalinity, dic, Ks)
    % [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = Ks.unpack();
    [K1,K2,KB] = Ks.unpack_some(["k1","k2","kb"]);
    boron = Ks.controls.composition.boron;

    g0 = K1.*K2.*KB.*(1-(2.*dic+boron)./alkalinity);
    g1 = K1.*(KB.*(1-boron./alkalinity-dic./alkalinity)+K2.*(1-2.*dic./alkalinity));
    g2 = KB.*(1-boron./alkalinity)+K1.*(1-dic./alkalinity);

    % Determine g21min
    g21min = g2.^2-3.*g1;
    g21min_positive = g21min > 0;
    sq21 = NaN(size(alkalinity,1),1);
    sq21(g21min_positive) = sqrt(g21min(g21min_positive));
    sq21(~g21min_positive) = 0;

    % Determine Hmin
    Hmin = NaN(size(alkalinity,1),1);
    g2_positive = g2 >=0;
    Hmin(~g2_positive) = (-g2(~g2_positive) + sq21(~g2_positive))./3;
    Hmin(g2_positive) = -g1(g2_positive)./(g2(g2_positive) + sq21(g2_positive));

    % Calculate initial pH
    pH_initial_guess = NaN(size(alkalinity,1),1);
    negative_alkalinity = alkalinity <= 0;
    pH_initial_guess(negative_alkalinity) = -log10(1e-3);

    medium_alkalinity = alkalinity > 0 & alkalinity < 2*dic + boron;
    pH_initial_guess(medium_alkalinity & g21min_positive) = ...
        -log10(Hmin(medium_alkalinity & g21min_positive) + ...
        sqrt(-(Hmin(medium_alkalinity & g21min_positive).^3 + g2(medium_alkalinity & g21min_positive).*Hmin(medium_alkalinity & g21min_positive).^2 + ...
        g1(medium_alkalinity & g21min_positive).*Hmin(medium_alkalinity & g21min_positive) + ...
        g0(medium_alkalinity & g21min_positive))./sq21(medium_alkalinity & g21min_positive)));
    pH_initial_guess(medium_alkalinity & ~g21min_positive) = -log10(1e-7);

    high_alkalinity = alkalinity >= 2.*dic + boron;
    pH_initial_guess(high_alkalinity) = -log10(1e-10);

    pH_out = pH_initial_guess;
end

function pH_out = calculate_pH_from_alkalinity_co2_munhoven(alkalinity, co2, Ks)
    [K1,K2,KB] = Ks.unpack_some(["K1","K2","KB"]);
    boron = Ks.controls.composition.boron;

    g0 = -2.*K1.*K2.*KB.*co2./alkalinity;
    g1 = -K1.*(2.*K2.*co2+KB.*co2)./alkalinity;
    g2 = KB-(boron.*KB+K1.*co2)./alkalinity;

    % Determine Hmin
    g21min = g2.^2-3.*g1;
    g21min_positive = g21min > 0;

    sq21 = NaN(size(alkalinity,1),1);
    sq21(g21min_positive) = sqrt(g21min(g21min_positive));
    sq21(~g21min_positive) = 0;

    g2_positive = g2 >=0;
    Hmin = NaN(size(alkalinity,1),1);
    Hmin(~g2_positive) = (-g2(~g2_positive) + sq21(~g2_positive))./3;
    Hmin(g2_positive) = -g1(g2_positive)./(g2(g2_positive) + sq21(g2_positive));

    % Calculate initial pH
    pH_initial_guess = NaN(size(alkalinity,1),1);

    negative_alkalinity = alkalinity <= 0;
    pH_initial_guess(negative_alkalinity) = -log10(1e-3);

    positive_alkalinity = alkalinity > 0;
    pH_initial_guess(positive_alkalinity & g21min_positive) = ...
        -log10(Hmin(positive_alkalinity & g21min_positive) + ...
        sqrt(-(Hmin(positive_alkalinity & g21min_positive).^3 + g2(positive_alkalinity & g21min_positive).*Hmin(positive_alkalinity & g21min_positive).^2 + ...
        g1(positive_alkalinity & g21min_positive).*Hmin(positive_alkalinity & g21min_positive)+...
        g0(positive_alkalinity & g21min_positive))./sq21(positive_alkalinity & g21min_positive)));

    pH_initial_guess(positive_alkalinity & ~g21min_positive) = -log10(1e-7);

    pH_out = pH_initial_guess;
end

function pH_out = calculate_pH_from_alkalinity_hco3_munhoven(alkalinity, hco3, Ks)
    [K2,KB] = Ks.unpack_some(["K2","KB"]);
    boron = Ks.controls.composition.boron;

    g0 = 2.*K2.*KB.*hco3;
    g1 = KB.*(hco3+boron-alkalinity)+2.*K2.*hco3;
    g2 = hco3-alkalinity;

    % Calculate initial pH
    pH_initial_guess = NaN(size(alkalinity,1),1);

    low_alkalinity = alkalinity <= hco3;
    pH_initial_guess(low_alkalinity) = -log10(1e-3);

    high_alkalinity = alkalinity > hco3;
    pH_initial_guess(high_alkalinity) = ...
        -log10((-g1(high_alkalinity)-sqrt(g1(high_alkalinity).^2-4.*g0(high_alkalinity).*g2(high_alkalinity)))./(2.*g2(high_alkalinity)));
    
    pH_out = pH_initial_guess;
end

function pH_out = calculate_pH_from_alkalinity_co3_munhoven(alkalinity, co3, Ks)
    [K2,KB] = Ks.unpack_some(["K2","KB"]);
    boron = Ks.controls.composition.boron;

    g0 = K2.*KB.*(2.*co3+boron-alkalinity);
    g1 = KB.*co3+K2.*(2.*co3-alkalinity);
    g2 = co3;

    % Calculate initial pH
    pH_initial_guess = NaN(size(alkalinity,1),1);

    low_alkalinity = alkalinity <= 2.*co3+boron;
    pH_initial_guess(low_alkalinity) = -log10(1e-3);

    high_alkalinity = alkalinity > 2.*co3+boron;
    pH_initial_guess(high_alkalinity) = ...
        -log10((-g1(high_alkalinity)+sqrt(g1(high_alkalinity).^2-4.*g0(high_alkalinity).*g2(high_alkalinity)))./(2.*g2(high_alkalinity)));
    
    pH_out = pH_initial_guess;
end


%% Revelle - buffering
function revelle = calculate_revelle_factor(alkalinity, dic,Ks)
    dic_start = dic;
    delta_dic = 0.00000001;% ' 0.01 umol/kg-SW (lower than prior versions of CO2SYS)


    % ' Find fCO2 at TA, TC + dTC
    dic_high = dic_start + delta_dic;
    pH_high = calculate_pH_from_alkalinity_dic(alkalinity, dic_high,Ks);
    fco2_high = calculate_fco2_from_dic_pH(dic_high, pH_high, Ks);

    % ' Find fCO2 at TA, TC - dTC
    dic_low = dic_start - delta_dic;
    pH_low = calculate_pH_from_alkalinity_dic(alkalinity, dic_low,Ks);
    fco2_low = calculate_fco2_from_dic_pH(dic_low, pH_low, Ks);

    % CalculateRevelleFactor:
    revelle = (fco2_high - fco2_low)./delta_dic./((fco2_high + fco2_low)./dic_start); % Corrected error pointed out by MP Humphreys (https://pyco2sys.readthedocs.io/en/latest/validate/)
end


%% Solubility
function [saturation_state_calcite,saturation_state_aragonite] = calculate_carbonate_solubility(salinity, temp_c, dic, pH, Ks, calcium,which_k1_k2,pressure)
    [K1,K2] = Ks.unpack_some(["K1","K2"]);

    temp_k = temp_c + 273.15;
    log_temp_k = log(temp_k);
    sqrt_salinity = sqrt(salinity);
    gas_constant = Constants.gas_constant;

    KAr = NaN(numel(temp_c),1);

    new_selected = (which_k1_k2~=6 & which_k1_k2~=7);
    if any(new_selected)
    % (below here, selected isn't used, since almost always all rows match the above criterium,
    %  in all other cases the rows will be overwritten later on).
        
        % AragoniteSolubility:
        % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKAr = (-171.945 - 0.077993.*temp_k(new_selected) + 2903.293./temp_k(new_selected)...
                  + 71.595.*log_temp_k(new_selected)./log(10)...
                  + (-0.068393 + 0.0017276.*temp_k(new_selected) + 88.135./temp_k(new_selected)).*sqrt_salinity(new_selected)...
                  - 0.10018.*salinity(new_selected) + 0.0059415.*sqrt_salinity(new_selected).*salinity(new_selected));
        KAr(new_selected)    = 10.^(logKAr);% ' this is in (mol/kg-SW)^2
                
        % PressureCorrectionForAragonite:
        % '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
        % '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
        % '       and 10^3 for Kappa factor)
        deltaVKAr = -48.76 + 0.5304.*temp_c(new_selected) + 2.8;
        KappaKAr  = (-11.76 + 0.3692.*temp_c(new_selected))./1000;
        lnKArfac  = (-deltaVKAr + 0.5.*KappaKAr.*pressure(new_selected)).*pressure(new_selected)./(gas_constant.*temp_k(new_selected));
        KAr(new_selected)       = KAr(new_selected).*exp(lnKArfac);
    end

    new_selected=(which_k1_k2==6 | which_k1_k2==7);
    if any(new_selected)
        % this is in (mol/kg-SW)^2
        %
        % *** CalculateKArforGEOSECS:
        % Berner, R. A., American Journal of Science 276:713-730, 1976:
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
        KAr(new_selected) = 1.45.*0.0000001.*(-34.452 - 39.866.*salinity(new_selected).^(1./3) +...
            110.21.*log(salinity(new_selected))./log(10) - 0.0000075752.*temp_k(new_selected).^2);% ' this is in (mol/kg-SW)^2
        % Berner (p. 722) states that he uses 1.48.
        % It appears that 1.45 was used in the GEOSECS calculations
        %
        % *** CalculatePressureEffectsOnKCaKArGEOSECS:
        % Culberson and Pytkowicz, Limnology and Oceanography 13:403-417, 1968
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
        % but their paper is not even on this topic).
        % The fits appears to be new in the GEOSECS report.
        % I can't find them anywhere else.
        KAr(new_selected) = KAr(new_selected).*exp((33.3 - 0.22.*temp_c(new_selected)).*pressure(new_selected)./(gas_constant.*temp_k(new_selected)));
    end
    
    % CalculateOmegasHere:
    H = 10.^(-pH);
    CO3 = dic.*K1.*K2./(K1.*H + H.*H + K1.*K2);
    saturation_state_calcite = CO3.*calcium./Ks.kc; % OmegaCa, dimensionless
    saturation_state_aragonite = CO3.*calcium./KAr; % OmegaAr, dimensionless
end
    

%% Utility
function [boron_alkalinity,hydroxide_alkalinity,phosphate_alkalinity,silicate_alkalinity,ammonia_alkalinity,sulphide_alkalinity,hydrogen_free,sulphate_acidity,fluorine_acidity] = calculate_alkalinity_parts(pH,Ks)
    [KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = Ks.unpack_some(["KW","KB","KF","KS","KP1","KP2","KP3","KSi","KNH4","KH2S"]);
    
    [ammonia,boron,fluorine,phosphate,silicate,sulphate,sulphide] = Ks.controls.composition.unpack_alkalinity();
    
    H = 10.^(-pH);
    boron_alkalinity = boron.*KB./(KB + H);
    hydroxide_alkalinity = KW./H;
    phosphate_alkalinity = phosphate.*(KP1.*KP2.*H + 2.*KP1.*KP2.*KP3 - H.*H.*H)./(H.*H.*H + KP1.*H.*H + KP1.*KP2.*H + KP1.*KP2.*KP3);
    silicate_alkalinity = silicate.*KSi./(KSi + H);
    ammonia_alkalinity = ammonia.*KNH4./(KNH4 + H);
    sulphide_alkalinity = sulphide.*KH2S./(KH2S + H);

    [~,~,pHfree,~] = Ks.controls.pH_scale_conversion(2).find_pH_on_all_scales(pH,Ks.controls);
        
    hydrogen_free = 10.^-pHfree; % this converts pHfree to Hfree
    sulphate_acidity = sulphate./(1 + KS./hydrogen_free); %' since KS is on the free scale
    fluorine_acidity = fluorine./(1 + KF./hydrogen_free); %' since KF is on the free scale
end

function show_failed(failed)
    if numel(failed)>1
        disp(['pH values did not converge for data on rows: ',num2str((failed)')]);
    else
        disp(['pH value did not converge for data on row: ',num2str((failed)')]);
    end
end
