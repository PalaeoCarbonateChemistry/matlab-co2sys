function [data,headers,nice_headers]=CO2SYS(parameter_1,parameter_2, ...
                                            parameter_1_type,parameter_2_type, ...
                                            salinity_in, ...
                                            temperature_in,temperature_out, ...
                                            pressure_in,pressure_out, ...
                                            silicate,phosphate,ammonia,sulphide, ...
                                            pH_scale_in, ...
                                            which_k1_k2,which_kso4,which_kf, which_boron, ...
                                            varargin)
    
    % Added by JM Epitalon
    % For computing derivative with respect to Ks, one has to call CO2sys with a perturbed K
    % Requested perturbation is passed through the following global variables
    global k_perturbation_GLOBAL    % Id of perturbed K
    global Perturb  % perturbation
    
    % Input conditioning
    
    % set default for optional input argument
    p_opt = 0;
    % parse optional input argument
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'co2_press')
            p_opt = varargin{i+1};
        end
    end
    
    % Determine lengths of input vectors
    veclengths=[length(parameter_1) length(parameter_2) length(parameter_1_type)...
                length(parameter_2_type) length(salinity_in) length(temperature_in)...
                length(temperature_out) length(pressure_in) length(pressure_out)...
                length(silicate) length(phosphate) length(ammonia) length(sulphide)...
                length(pH_scale_in) length(which_k1_k2) length(which_kso4)...
	            length(which_kf) length(which_boron)];
    
    if length(unique(veclengths))>2
	    disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end
    
    % Find the longest column vector:
    number_of_points = max(veclengths);
    
    % Populate column vectors
    parameter_1(1:number_of_points,1)          = parameter_1(:)          ;
    parameter_2(1:number_of_points,1)          = parameter_2(:)          ;
    parameter_1_type(1:number_of_points,1)      = parameter_1_type(:)      ;
    parameter_2_type(1:number_of_points,1)      = parameter_2_type(:)      ;
    salinity_in(1:number_of_points,1)           = salinity_in(:)           ;
    temperature_in(1:number_of_points,1)        = temperature_in(:)        ;
    temperature_out(1:number_of_points,1)       = temperature_out(:)       ;
    pressure_in(1:number_of_points,1)        = pressure_in(:)        ;
    pressure_out(1:number_of_points,1)       = pressure_out(:)       ;
    silicate(1:number_of_points,1)            = silicate(:)            ;
    phosphate(1:number_of_points,1)           = phosphate(:)           ;
    ammonia(1:number_of_points,1)           = ammonia(:)           ;
    sulphide(1:number_of_points,1)           = sulphide(:)           ;
    pH_scale_in(1:number_of_points,1)     = pH_scale_in(:)     ;
    which_k1_k2(1:number_of_points,1) = which_k1_k2(:) ;
    which_kso4(1:number_of_points,1)  = which_kso4(:)  ;
    which_kf(1:number_of_points,1)    = which_kf(:)    ;
    which_boron(1:number_of_points,1)         = which_boron(:)         ;
    
    % Assign input to the 'historical' variable names.
    % pH_scale_in_GLOBAL      = pH_scale_in;
    
    gas_constant = 83.14462618; % ml bar-1 K-1 mol-1,
    
    % Generate empty vectors for...
    TA   = nan(number_of_points,1); % Talk
    TC   = nan(number_of_points,1); % DIC
    PH   = nan(number_of_points,1); % pH
    PC   = nan(number_of_points,1); % pCO2
    FC   = nan(number_of_points,1); % fCO2
    HCO3 = nan(number_of_points,1); % [HCO3]
    CO3  = nan(number_of_points,1); % [CO3]
    CO2  = nan(number_of_points,1); % [CO2*]

    salinity = salinity_in;
    
    % Assign values to empty vectors.
    selected=(parameter_1_type==1 & parameter_1~=-999);   
    TA(selected)=parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg

    selected=(parameter_1_type==2 & parameter_1~=-999);   
    TC(selected)=parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg

    selected=(parameter_1_type==3 & parameter_1~=-999);   
    PH(selected)=parameter_1(selected);

    selected=(parameter_1_type==4 & parameter_1~=-999);   
    PC(selected)=parameter_1(selected)/1e6; % Convert from microatm. to atm.

    selected=(parameter_1_type==5 & parameter_1~=-999);   
    FC(selected)=parameter_1(selected)/1e6; % Convert from microatm. to atm.

    selected=(parameter_1_type==6 & parameter_1~=-999); 
    HCO3(selected)=parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_1_type==7 & parameter_1~=-999);  
    CO3(selected)=parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_1_type==8 & parameter_1~=-999);  
    CO2(selected)=parameter_1(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_2_type==1 & parameter_2~=-999);   
    TA(selected)=parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_2_type==2 & parameter_2~=-999);   
    TC(selected)=parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_2_type==3 & parameter_2~=-999);   
    PH(selected)=parameter_2(selected);
    
    selected=(parameter_2_type==4 & parameter_2~=-999);   
    PC(selected)=parameter_2(selected)/1e6; % Convert from microatm. to atm.
    
    selected=(parameter_2_type==5 & parameter_2~=-999);   
    FC(selected)=parameter_2(selected)/1e6; % Convert from microatm. to atm.
    
    selected=(parameter_2_type==6 & parameter_2~=-999); 
    HCO3(selected)=parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_2_type==7 & parameter_2~=-999);  
    CO3(selected)=parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    selected=(parameter_2_type==8 & parameter_2~=-999);  
    CO2(selected)=parameter_2(selected)/1e6; % Convert from micromol/kg to mol/kg
    
    % Generate the columns holding Si, Phos, Amm, H2S and salinity.
    % Pure Water case:
    selected=(which_k1_k2==8);
    salinity(selected) = 0;

    boron_concentration = calculate_boron_concentration(salinity,number_of_points,which_boron,which_k1_k2);
    fluorine_concentration = calculate_fluorine_concentration(salinity);
    sulphate_concentration = calculate_sulphate_concentration(salinity);
    calcium_concentration = calculate_calcium_concentration(salinity,number_of_points,which_k1_k2);
    phosphate_concentration = calculate_phosphate_concentration(phosphate,number_of_points,which_k1_k2);
    silicate_concentration = calculate_silicate_concentration(silicate,number_of_points,which_k1_k2);
    ammonia_concentration = calculate_ammonia_concentration(ammonia,number_of_points,which_k1_k2);
    sulphide_concentration = calculate_sulphide_concentration(sulphide,number_of_points,which_k1_k2);
  
    % The vector 'peng_correction' is used to modify the value of TA, for those
    % cases where which_k1_k2==7, since PAlk(Peng) = PAlk(Dickson) + phosphate.
    % Thus, peng_correction is 0 for all cases where which_k1_k2 is not 7
    peng_correction=zeros(number_of_points,1); 
    selected=which_k1_k2==7; 
    peng_correction(selected)=phosphate_concentration(selected);
    
    % Calculate the constants for all samples at input conditions
    % The constants calculated for each sample will be on the appropriate pH scale!
    Ks_in = calculate_equilibrium_constants(number_of_points,temperature_in,pressure_in/10,salinity,pH_scale_in,p_opt,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2);
    
    % Added by JM Epitalon
    % For computing derivative with respect to Ks, one has to perturb the value of one K
    % Requested perturbation is passed through global variables k_perturbation_GLOBAL and Perturb
    if (~ isempty(k_perturbation_GLOBAL))
        switch k_perturbation_GLOBAL
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
                boron_concentration = boron_concentration + Perturb;
        end
    end
    
    K0_in = Ks_in("K0");
    
    % Make sure fCO2 is available for each sample that has pCO2 or CO2.
    temp_k = temperature_in+273.15;
    fugacity_factor = calculate_fugacity_factor(p_opt,gas_constant,number_of_points,which_k1_k2,temp_k);
    
    selected = (~isnan(PC) & (parameter_1_type==4 | parameter_2_type==4));  
    FC(selected) = PC(selected).*fugacity_factor(selected);

    selected = (~isnan(CO2) & (parameter_1_type==8 | parameter_2_type==8)); 
    FC(selected) = CO2(selected)./K0_in(selected);
    
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
    combination = 10*min(parameter_1_type,parameter_2_type) + max(parameter_1_type,parameter_2_type);
    
    % Calculate missing values for AT,CT,PH,FC,HCO3,CO3,CO2:
    % pCO2 will be calculated later on, routines work with fCO2.
    selected=combination==12; % input TA, TC
    if any(selected)
    selected=(~isnan(TAc) & ~isnan(TCc) & selected);
        PHic(selected)                = CalculatepHfromTATC(TAc(selected)-peng_correction(selected),TCc(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
        selected=(~isnan(PHic) & selected);
        if any(selected)
           FCic(selected)              = CalculatefCO2fromTCpH(TCc(selected), PHic(selected),Ks_in,selected);
           [CO3ic(selected),HCO3ic(selected)] = CalculateCO3HCO3fromTCpH(TCc(selected),PHic(selected),Ks_in,selected);
        end
    end
    selected=combination==13; % input TA, pH
    if any(selected)
    selected=(~isnan(TAc) & ~isnan(PHic) & selected);
        TCc(selected)                 = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
        FCic(selected)                = CalculatefCO2fromTCpH(TCc(selected),PHic(selected),Ks_in,selected);
        [CO3ic(selected),HCO3ic(selected)]   = CalculateCO3HCO3fromTCpH(TCc(selected),PHic(selected),Ks_in,selected);
    end
    selected=combination==14 | combination==15 | combination==18; % input TA, (pCO2 or fCO2 or CO2)
    if any(selected)
    selected=(~isnan(TAc) & ~isnan(FCic) & selected);
        PHic(selected)                = CalculatepHfromTAfCO2(TAc(selected)-peng_correction(selected),FCic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
        selected=(~isnan(PHic) & selected);
        if any(selected)
           TCc(selected)              = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
           [CO3ic(selected),HCO3ic(selected)]= CalculateCO3HCO3fromTCpH(TCc(selected),PHic(selected),Ks_in,selected);
        end
    end
    selected=combination==16; % input TA, HCO3
    if any(selected)
    selected=(~isnan(TAc) & ~isnan(HCO3ic) & selected);
        PHic(selected)                = CalculatepHfromTAHCO3(TAc(selected)-peng_correction(selected),HCO3ic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);  % added Peng correction // MPH
        selected=(~isnan(PHic) & selected);
        if any(selected)
           TCc(selected)              = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
           FCic(selected)             = CalculatefCO2fromTCpH(TCc(selected),PHic(selected), Ks_in,selected); 
           CO3ic(selected)            = CalculateCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
        end
    end
    selected=combination==17; % input TA, CO3
    if any(selected)
    selected=(~isnan(TAc) & ~isnan(CO3ic) & selected);
        PHic(selected)                 = CalculatepHfromTACO3(TAc(selected)-peng_correction(selected),CO3ic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);  % added Peng correction // MPH
        selected=(~isnan(PHic) & selected);
        if any(selected)
           TCc(selected)               = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
           FCic(selected)              = CalculatefCO2fromTCpH(TCc(selected),PHic(selected), Ks_in,selected); 
           HCO3ic(selected)            = CalculateHCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
        end
    end
    selected=combination==23; % input TC, pH
    if any(selected)
    selected=(~isnan(TCc) & ~isnan(PHic) & selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        FCic(selected)                 = CalculatefCO2fromTCpH(TCc(selected),PHic(selected),Ks_in,selected);
        [CO3ic(selected),HCO3ic(selected)]    = CalculateCO3HCO3fromTCpH(TCc(selected), PHic(selected),Ks_in,selected);
    end
    selected=combination==24 | combination==25 | combination==28;  % input TC, (pCO2 or fCO2 or CO2)
    if any(selected)
    selected=(~isnan(TCc) & ~isnan(FCic) & selected);
        PHic(selected)                 = CalculatepHfromTCfCO2(TCc(selected),FCic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        [CO3ic(selected),HCO3ic(selected)]    = CalculateCO3HCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==26; % input TC, HCO3
    if any(selected)
    selected=(~isnan(TCc) & ~isnan(HCO3ic) & selected);
        [PHic(selected),FCic(selected)]       = CalculatepHfCO2fromTCHCO3(TCc(selected),HCO3ic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        CO3ic(selected)                = CalculateCO3fromTCpH(TCc(selected),PHic(selected),Ks_in,selected);
    end
    selected=combination==27; % input TC, CO3
    if any(selected)
    selected=(~isnan(TCc) & ~isnan(CO3ic) & selected);
        [PHic(selected),FCic(selected)]       = CalculatepHfCO2fromTCCO3(TCc(selected),CO3ic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        HCO3ic(selected)               = CalculateHCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==34 | combination==35 | combination==38; % input pH, (pCO2 or fCO2 or CO2)
    if any(selected)
    selected=(~isnan(PHic) & ~isnan(FCic) & selected);
        TCc(selected)                  = CalculateTCfrompHfCO2(PHic(selected),FCic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        [CO3ic(selected),HCO3ic(selected)]    = CalculateCO3HCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==36; % input pH, HCO3
    if any(selected)
    selected=(~isnan(PHic) & ~isnan(HCO3ic) & selected);
        TAc(selected)                  = CalculateTAfrompHHCO3(PHic(selected),HCO3ic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        TCc(selected)                  = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
        FCic(selected)                 = CalculatefCO2fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
        CO3ic(selected)                = CalculateCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==37; % input pH, CO3
    if any(selected)
    selected=(~isnan(PHic) & ~isnan(CO3ic) & selected);
        TAc(selected)                  = CalculateTAfrompHCO3(PHic(selected),CO3ic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        TCc(selected)                  = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
        FCic(selected)                 = CalculatefCO2fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
        HCO3ic(selected)               = CalculateHCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==46 | combination==56 | combination==68; % input (pCO2 or fCO2 or CO2), HCO3
    if any(selected)
    selected=(~isnan(FCic) & ~isnan(HCO3ic) & selected);
        PHic(selected)                 = CalculatepHfromfCO2HCO3(FCic(selected),HCO3ic(selected), Ks_in,selected);
        TCc(selected)                  = CalculateTCfrompHfCO2(PHic(selected),FCic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        CO3ic(selected)                = CalculateCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==47 | combination==57 | combination==78; % input (pCO2 or fCO2 or CO2), CO3
    if any(selected)
    selected=(~isnan(FCic) & ~isnan(CO3ic) & selected);
        PHic(selected)                 = CalculatepHfromfCO2CO3(FCic(selected),CO3ic(selected), Ks_in,selected);
        TCc(selected)                  = CalculateTCfrompHfCO2 (PHic(selected),FCic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfromTCpH(TCc(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        HCO3ic(selected)               = CalculateHCO3fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
    end
    selected=combination==67; % input HCO3, CO3
    if any(selected)
    selected=(~isnan(HCO3ic) & ~isnan(CO3ic) & selected);
        PHic(selected)                 = CalculatepHfromCO3HCO3(CO3ic(selected),HCO3ic(selected), Ks_in,selected);
        TAc(selected)                  = CalculateTAfrompHCO3(PHic(selected),CO3ic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k) + peng_correction(selected);
        TCc(selected)                  = CalculateTCfromTApH(TAc(selected)-peng_correction(selected),PHic(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
        FCic(selected)                 = CalculatefCO2fromTCpH(TCc(selected),PHic(selected), Ks_in,selected);
        %CO2ic(selected)                = CalculateCO2fromTCpH(TCc(selected),PHic(selected));
    end
    
    % By now, an fCO2 value is available for each sample.
    % Generate the associated pCO2 value:
    selected = (isnan(PCic) & (parameter_1_type~=4 | parameter_2_type~=4)); PCic(selected)  = FCic(selected)./fugacity_factor(selected);
    % Generate the associated CO2 value:
    selected = (isnan(CO2ic) & (parameter_1_type~=8 | parameter_2_type~=8)); CO2ic(selected) = FCic(selected).*K0_in(selected);
    
    % Calculate Other Params At Input Conditions:
    BAlkinp    = nan(number_of_points,1); % Generate empty vectors
    [OHinp,PAlkinp,SiAlkinp,AmmAlkinp,HSAlkinp,Hfreeinp,HSO4inp,HFinp,...
        Revelleinp,OmegaCainp,OmegaArinp,xCO2dryinp] = deal(BAlkinp);
    selected=(~isnan(PHic)); % if PHic = NaN, pH calculation was not performed or did not converge
    [BAlkinp(selected),OHinp(selected), PAlkinp(selected),SiAlkinp(selected),AmmAlkinp(selected),...
        HSAlkinp(selected), Hfreeinp(selected),HSO4inp(selected),HFinp(selected)] = CalculateAlkParts(PHic,pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    PAlkinp(selected)                = PAlkinp(selected)+peng_correction(selected);
    Revelleinp(selected)             = RevelleFactor(TAc(selected)-peng_correction(selected), TCc(selected),pH_scale_in,Ks_in,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    [OmegaCainp(selected),OmegaArinp(selected)] = CaSolubility(salinity(selected), temperature_in(selected), TCc(selected), PHic(selected), Ks_in, sqrt(salinity(selected)),gas_constant,calcium_concentration,which_k1_k2,pressure_in/10,selected);
    VPFac = calculate_VPFac(salinity,temp_k);
    xCO2dryinp(~isnan(PCic),1) = PCic(~isnan(PCic),1)./VPFac(~isnan(PCic),1); % ' this assumes pTot = 1 atm
    SIRinp = HCO3ic./(Hfreeinp.*1e6);
    
    % % Just for reference, convert pH at input conditions to the other scales
    pHicT = nan(number_of_points,1);
    pHicS = nan(number_of_points,1);
    pHicF = nan(number_of_points,1);
    pHicN = nan(number_of_points,1);
    [pHicT(selected),pHicS(selected),pHicF(selected),pHicN(selected)]=FindpHOnAllScales(PHic(selected),pH_scale_in,Ks_in,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    
    % Merge the Ks at input into an array. Ks at output will be glued to this later.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks_in);
    KIVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S];

    clear K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S
    
    % Calculate the constants for all samples at output conditions
    Ks_out = calculate_equilibrium_constants(number_of_points,temperature_out,pressure_out/10,salinity,pH_scale_in,p_opt,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2);
    
    % Added by JM Epitalon
    % For computing derivative with respect to Ks, one has to perturb the value of one K
    % Requested perturbation is passed through global variables k_perturbation_GLOBAL and Perturb
    if (~ isempty(k_perturbation_GLOBAL))
        switch k_perturbation_GLOBAL
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
                boron_concentration = boron_concentration + Perturb;
        end
    end                  

    % For output conditions, using conservative TA and TC, calculate pH, fCO2
    % and pCO2, HCO3, CO3, and CO2
    temp_k = temperature_out+273.15;
    selected=(~isnan(TAc) & ~isnan(TCc)); % i.e., do for all samples that have TA and TC values
    PHoc=nan(number_of_points,1);
    [CO3oc,HCO3oc,FCoc] = deal(PHoc);
    PHoc(selected) = CalculatepHfromTATC(TAc(selected)-peng_correction(selected), TCc(selected),pH_scale_in,Ks_out,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
        FCoc(selected) = CalculatefCO2fromTCpH(TCc(selected), PHoc(selected), Ks_out,selected);
        [CO3oc(selected),HCO3oc(selected)] = CalculateCO3HCO3fromTCpH(TCc(selected),PHoc(selected), Ks_out,selected);
    
    % Generate the associated pCO2 value:
    fugacity_factor = calculate_fugacity_factor(p_opt,gas_constant,number_of_points,which_k1_k2,temp_k);
    PCoc  = FCoc./fugacity_factor;
    % Generate the associated CO2 value:

    K0_out = Ks_out("K0");
    CO2oc = FCoc.*K0_out;
    
    % Calculate Other Params At Output Conditions
    BAlkout    = nan(number_of_points,1); % Generate empty vectors
    [OHout,PAlkout,SiAlkout,AmmAlkout,HSAlkout,Hfreeout,HSO4out,HFout,...
        Revelleout,OmegaCaout,OmegaArout,xCO2dryout] = deal(BAlkout);
    selected=(~isnan(PHoc)); % if PHoc = NaN, pH calculation was not performed or did not converge
    [BAlkout(selected),OHout(selected),PAlkout(selected),SiAlkout(selected),AmmAlkout(selected),...
        HSAlkout(selected), Hfreeout(selected),HSO4out(selected),HFout(selected)] = CalculateAlkParts(PHoc,pH_scale_in,Ks_out,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    PAlkout(selected)                 = PAlkout(selected)+peng_correction(selected);
    Revelleout(selected)              = RevelleFactor(TAc(selected)-peng_correction(selected), TCc(selected),pH_scale_in,Ks_out,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    [OmegaCaout(selected),OmegaArout(selected)] = CaSolubility(salinity(selected), temperature_out(selected), TCc(selected), PHoc(selected), Ks_out, sqrt(salinity(selected)), gas_constant, calcium_concentration,which_k1_k2,pressure_out/10,selected);
    VPFac = calculate_VPFac(salinity,temp_k);
    xCO2dryout(~isnan(PCoc),1)    = PCoc(~isnan(PCoc))./VPFac(~isnan(PCoc)); % ' this assumes pTot = 1 atm
    SIRout = HCO3oc./(Hfreeout.*1e6);
    
    % Just for reference, convert pH at output conditions to the other scales
    pHocT = nan(number_of_points,1);
    pHocS = nan(number_of_points,1);
    pHocF = nan(number_of_points,1);
    pHocN = nan(number_of_points,1);
    [pHocT(selected),pHocS(selected),pHocF(selected),pHocN(selected)]=FindpHOnAllScales(PHoc(selected),pH_scale_in,Ks_out,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks_out);
    KOVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S];
    TVEC =[boron_concentration fluorine_concentration sulphate_concentration phosphate_concentration silicate_concentration ammonia_concentration sulphide_concentration];
    
    % Saving data in array, 99 columns, as many rows as samples input
    data=[TAc*1e6         TCc*1e6        PHic           PCic*1e6        FCic*1e6...
          HCO3ic*1e6      CO3ic*1e6      CO2ic*1e6      BAlkinp*1e6     OHinp*1e6...
          PAlkinp*1e6     SiAlkinp*1e6   AmmAlkinp*1e6  HSAlkinp*1e6    Hfreeinp*1e6... %%% Multiplied Hfreeinp *1e6, svh20100827
          Revelleinp      OmegaCainp     OmegaArinp     xCO2dryinp*1e6  SIRinp...
          PHoc            PCoc*1e6       FCoc*1e6       HCO3oc*1e6      CO3oc*1e6...
          CO2oc*1e6       BAlkout*1e6    OHout*1e6      PAlkout*1e6     SiAlkout*1e6...
          AmmAlkout*1e6   HSAlkout*1e6   Hfreeout*1e6   Revelleout      OmegaCaout... %%% Multiplied Hfreeout *1e6, svh20100827
          OmegaArout      xCO2dryout*1e6 SIRout         pHicT           pHicS...
          pHicF           pHicN          pHocT          pHocS           pHocF...
          pHocN           temperature_in         temperature_out        pressure_in          pressure_out...
          parameter_1_type        parameter_2_type       which_k1_k2  which_kso4    which_kf...
          which_boron           pH_scale_in      salinity_in            phosphate             silicate...
          ammonia             sulphide            KIVEC          KOVEC           TVEC*1e6];
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
    
    clear global K2 KP3 
    clear global KB KS T BORON 
    clear global K KF KSi KNH4 KH2S 
    clear global K0 KP1 KW 
    clear global K1 KP2
	
end 

function varargout=CalculatepHfromTATC(TAi, TCi, pH_scale_in, Ks, boron_concentration, fluorine_concentration, sulphate_concentration, phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
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

    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    K1F=K1(selected);     K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
    vl          = sum(selected);  % VectorLength
    % Find initital pH guess using method of Munhoven (2013)
    pHGuess         = CalculatepHfromTATCMunhoven(TAi, TCi, Ks, boron_concentration,selected);
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
        [~,~,pHfree,~] = FindpHOnAllScales(pH,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
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
end

function varargout=CalculatefCO2fromTCpH(TCx, pHx, Ks,selected)
    % ' SUB CalculatefCO2fromTCpH, version 02.02, 12-13-96, written by Ernie Lewis.
    % ' Inputs: TC, pH, K0, K1, K2
    % ' Output: fCO2
    % ' This calculates fCO2 from TC and pH, using K0, K1, and K2.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    H            = 10.^(-pHx);
    fCO2x        = TCx.*H.*H./(H.*H + K1(selected).*H + K1(selected).*K2(selected))./K0(selected);
    varargout{1} = fCO2x;
    end % end nested function
    
function varargout=CalculateTCfromTApH(TAx, pHx,pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    
    K1F=K1(selected);     K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
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
    [~,~,pHfree,~] = FindpHOnAllScales(pHx,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
    Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree); %' since KF is on the free scale
    CAlk      = TAx - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    TCxtemp   = CAlk.*(H.*H + K1F.*H + K1F.*K2F)./(K1F.*(H + 2.*K2F));
    varargout{1} = TCxtemp;
end % end nested function

function varargout=CalculatepHfromTAfCO2(TAi, fCO2i, pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
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
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    K0F=K0(selected);     K1F=K1(selected);     K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
    vl         = sum(selected); % vectorlength
    % Find initital pH guess using method of Munhoven (2013)
    CO2i       = fCO2i.*K0F; % Convert fCO2 to CO2
    pHGuess    = CalculatepHfromTACO2Munhoven(TAi, CO2i, Ks,boron_concentration,selected);
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
        [~,~,pHfree,~] = FindpHOnAllScales(pH,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
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
end

function varargout=CalculateTAfromTCpH(TCi, pHi,pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    % ' SUB CalculateTAfromTCpH, version 02.02, 10-10-97, written by Ernie Lewis.
    % ' Inputs: TC, pH, K(), T()
    % ' Output: TA
    % ' This calculates TA from TC and pH.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    K1F=K1(selected);     K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
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
    [~,~,pHfree,~] = FindpHOnAllScales(pHi,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    TActemp    = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
    varargout{1}=TActemp;
end

function varargout=CalculatepHfromTCfCO2(TCi, fCO2i, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculatepHfromTCfCO2, version 02.02, 11-12-96, written by Ernie Lewis.
    % ' Inputs: TC, fCO2, K0, K1, K2
    % ' Output: pH
    % ' This calculates pH from TC and fCO2 using K0, K1, and K2 by solving the
    % '       quadratic in H: fCO2.*K0 = TC.*H.*H./(K1.*H + H.*H + K1.*K2).
    % ' if there is not a real root, then pH is returned as missingn.
    RR = K0(selected).*fCO2i./TCi;
    %       if RR >= 1
    %          varargout{1}= missingn;
    %          disp('nein!');return;
    %       end
    % check after sub to see if pH = missingn.
    Discr = (K1(selected).*RR).*(K1(selected).*RR) + 4.*(1 - RR).*(K1(selected).*K2(selected).*RR);
    H     = 0.5.*(K1(selected).*RR + sqrt(Discr))./(1 - RR);
    %       if (H <= 0)
    %           pHctemp = missingn;
    %       else
    pHctemp = log(H)./log(0.1);
    %       end
    varargout{1}=pHctemp;
end

function varargout=CalculateTCfrompHfCO2(pHi, fCO2i, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculateTCfrompHfCO2, version 01.02, 12-13-96, written by Ernie Lewis.
    % ' Inputs: pH, fCO2, K0, K1, K2
    % ' Output: TC
    % ' This calculates TC from pH and fCO2, using K0, K1, and K2.
    H       = 10.^(-pHi);
    TCctemp = K0(selected).*fCO2i.*(H.*H + K1(selected).*H + K1(selected).*K2(selected))./(H.*H);
    varargout{1}=TCctemp;
end

function varargout=CalculateTAfrompHHCO3(pHi, HCO3i,pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    % ' SUB CalculateTAfrompHCO3, version 01.0, 3-19, added by J. Sharp
    % ' Inputs: pH, HCO3, K(), T()
    % ' Output: TA
    % ' This calculates TA from pH and HCO3.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    K1F=K1(selected);     K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
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
    [~,~,pHfree,~] = FindpHOnAllScales(pHi,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    TActemp     = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
    varargout{1}=TActemp;
end

function varargout=CalculatepHfromTAHCO3(TAi, HCO3i,pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
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
    K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
    vl         = sum(selected); % vectorlength
    % Find initital pH guess using method of Munhoven (2013)
    pHGuess    = CalculatepHfromTAHCO3Munhoven(TAi, HCO3i, Ks,boron_concentration,selected);
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
        [~,~,pHfree,~] = FindpHOnAllScales(pH,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
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
end

function varargout=CalculatepHfromTCHCO3(TCi, HCO3i, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
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
    Discr = ((1-RR).*(1-RR) - 4.*(1./(K1(selected))).*(K2(selected)));
    H     = 0.5.*((-(1-RR)) - sqrt(Discr))./(1./(K1(selected))); % Subtraction
    %       if (H <= 0)
    %           pHctemp = missingn;
    %       else
    pHctemp = log(H)./log(0.1);
    %       end
    varargout{1}=pHctemp;
end

function varargout=CalculatepHfromfCO2HCO3(fCO2i, HCO3i, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculatepHfromfCO2HCO3, version 01.0, 3-19, added by J. Sharp
    % ' Inputs: fCO2, HCO3, K0, K1, K2
    % ' Output: pH
    % ' This calculates pH from fCO2 and HCO3, using K0, K1, and K2.
    H            = (fCO2i.*K0(selected).*K1(selected))./HCO3i;  % removed incorrect (selected) index from HCO3i // MPH
    pHx          = -log10(H);
    varargout{1} = pHx;
    end % end nested function
    
    function varargout=CalculatepHfCO2fromTCHCO3(TCx, HCO3x, Ks,selected)
    % Outputs pH fCO2, in that order
    % SUB CalculatepHfCO2fromTCHCO3, version 01.0, 3-19, added by J. Sharp
    % Inputs: pHScale%, which_k1_k2%, which_kso4%, TC, HCO3, salinity, K(), T(), TempC, Pdbar
    % Outputs: pH, fCO2
    % This calculates pH and fCO2 from TC and HCO3 at output conditions.
    pHx   = CalculatepHfromTCHCO3(TCx, HCO3x, Ks,selected); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
    fCO2x = CalculatefCO2fromTCpH(TCx, pHx, Ks,selected);
    varargout{1} = pHx;
    varargout{2} = fCO2x;
    end


    function varargout=CalculateTAfrompHCO3(pHi, CO3i,pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculateTAfrompHCO3, version 01.0, 8-18, added by J. Sharp
    % ' Inputs: pH, CO3, K(), T()
    % ' Output: TA
    % ' This calculates TA from pH and CO3.
    K1F=K1(selected);     K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
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
    [~,~,pHfree,~] = FindpHOnAllScales(pHi,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    TActemp     = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
    varargout{1}=TActemp;
end

function varargout=CalculatepHfromTACO3(TAi,CO3i,pH_scale_in,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
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
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    K2F=K2(selected);     KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
    vl         = sum(selected); % vectorlength
    % Find initital pH guess using method of Munhoven (2013)
    pHGuess    = CalculatepHfromTACO3Munhoven(TAi, CO3i, Ks,boron_concentration,selected);
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
        [~,~,pHfree,~] = FindpHOnAllScales(pH,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
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
end

function varargout=CalculatepHfromTCCO3(TCi, CO3i, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
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
    Discr = ((1./K2(selected)).*(1./K2(selected)) - 4.*(1./(K1(selected).*K2(selected))).*(1-RR));
    H     = 0.5.*((-1./K2(selected)) + sqrt(Discr))./(1./(K1(selected).*K2(selected))); % Addition
    %       if (H <= 0)
    %           pHctemp = missingn;
    %       else
    pHctemp = log(H)./log(0.1);
    %       end
    varargout{1}=pHctemp;
end

function varargout=CalculatepHfromfCO2CO3(fCO2i, CO3i, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculatepHfromfCO2CO3, version 01.0, 8-18, added by J. Sharp
    % ' Inputs: fCO2, CO3, K0, K1, K2
    % ' Output: pH
    % ' This calculates pH from fCO2 and CO3, using K0, K1, and K2.
    H            = sqrt((fCO2i.*K0(selected).*K1(selected).*K2(selected))./CO3i);    % removed incorrect (selected) index from CO3i // MPH
    pHx          = -log10(H);
    varargout{1} = pHx;
end

function varargout=CalculatepHfCO2fromTCCO3(TCx, CO3x, Ks,selected)
    % Outputs pH fCO2, in that order
    % SUB CalculatepHfCO2fromTCCO3, version 01.0, 8-18, added by J. Sharp
    % Inputs: pHScale%, which_k1_k2%, which_kso4%, TC, CO3, salinity, K(), T(), TempC, Pdbar
    % Outputs: pH, fCO2
    % This calculates pH and fCO2 from TC and CO3 at output conditions.
    pHx   = CalculatepHfromTCCO3(TCx, CO3x, Ks,selected); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
    fCO2x = CalculatefCO2fromTCpH(TCx, pHx, Ks,selected);
    varargout{1} = pHx;
    varargout{2} = fCO2x;
end

function varargout=CalculatepHfromCO3HCO3(CO3x, HCO3x, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculatepHfromCO3HCO3, version 01.0, 3-19, added by J. Sharp
    % ' Inputs: CO3, HCO3, K2
    % ' Output: pH
    % ' This calculates fCO2 from TC and pH, using K2.
    H            = HCO3x.*K2(selected)./CO3x;
    pHx          = -log10(H);
    varargout{1} = pHx;
end

function varargout=CalculateCO3HCO3fromTCpH(TCx, pHx, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculateCO3HCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
    % ' Inputs: TC, pH, K1, K2
    % ' Output: CO3, HCO3, CO2
    % ' This calculates CO3, HCO3, and CO2 from TC and pH, using K1, and K2.
    H            = 10.^(-pHx);
    CO3x         = TCx.*K1(selected).*K2(selected)./(K1(selected).*H + H.*H + K1(selected).*K2(selected));
    HCO3x        = TCx.*K1(selected).*H./(K1(selected).*H + H.*H + K1(selected).*K2(selected));
    varargout{1} = CO3x;
    varargout{2} = HCO3x;
end

function varargout=CalculateCO3fromTCpH(TCx, pHx, Ks,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    % ' SUB CalculateCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
    % ' Inputs: TC, pH, K1, K2
    % ' Output: CO3, CO2
    % ' This calculates CO3 and CO2 from TC and pH, using K1, and K2.
    H            = 10.^(-pHx);
    CO3x         = TCx.*K1(selected).*K2(selected)./(K1(selected).*H + H.*H + K1(selected).*K2(selected));
    varargout{1} = CO3x;
end

function varargout=CalculateHCO3fromTCpH(TCx, pHx, Ks,selected)
    % ' SUB CalculateHCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
    % ' Inputs: TC, pH, K1, K2
    % ' Output: HCO3, CO2
    % ' This calculates HCO3 and CO2 from TC and pH, using K1, and K2.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    H            = 10.^(-pHx);
    HCO3x        = TCx.*K1(selected).*H./(K1(selected).*H + H.*H + K1(selected).*K2(selected));
    varargout{1} = HCO3x;
end


function varargout=CalculatepHfromTATCMunhoven(TAi, TCi, Ks,boron_concentration,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    K1F=K1(selected);     K2F=K2(selected);     TBF =boron_concentration(selected);    KBF=KB(selected);
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

function varargout=CalculatepHfromTACO2Munhoven(TAi, CO2x, Ks, boron_concentration,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    K1F=K1(selected);     K2F=K2(selected);     TBF =boron_concentration(selected);    KBF=KB(selected);
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

function varargout=CalculatepHfromTAHCO3Munhoven(TAi, HCO3x, Ks, boron_concentration,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    K1F=K1(selected);     K2F=K2(selected);     TBF =boron_concentration(selected);    KBF=KB(selected);
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

function varargout=CalculatepHfromTACO3Munhoven(TAi, CO3x, Ks, boron_concentration,selected)
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    K1F=K1(selected);     K2F=K2(selected);     TBF =boron_concentration(selected);    KBF=KB(selected);
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

function varargout=RevelleFactor(TAi, TCi, pH_scale,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    % ' SUB RevelleFactor, version 01.03, 01-07-97, written by Ernie Lewis.
    % ' Inputs: which_k1_k2%, TA, TC, K0, K(), T()
    % ' Outputs: Revelle
    % ' This calculates the Revelle factor (dfCO2/dTC)|TA/(fCO2/TC).
    % ' It only makes sense to talk about it at pTot = 1 atm, but it is computed
    % '       here at the given K(), which may be at pressure <> 1 atm. Care must
    % '       thus be used to see if there is any validity to the number computed.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    TC0 = TCi;
    dTC = 0.00000001;% ' 0.01 umol/kg-SW (lower than prior versions of CO2SYS)
    % ' Find fCO2 at TA, TC + dTC
    TCi = TC0 + dTC;
    pHc= CalculatepHfromTATC(TAi, TCi, pH_scale,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    fCO2c= CalculatefCO2fromTCpH(TCi, pHc, Ks,selected);
    fCO2plus = fCO2c;
    % ' Find fCO2 at TA, TC - dTC
    TCi = TC0 - dTC;
    pHc= CalculatepHfromTATC(TAi, TCi, pH_scale,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k);
    fCO2c= CalculatefCO2fromTCpH(TCi, pHc, Ks,selected);
    fCO2minus = fCO2c;
    % CalculateRevelleFactor:
    Revelle = (fCO2plus - fCO2minus)./dTC./((fCO2plus + fCO2minus)./TC0); % Corrected error pointed out by MP Humphreys (https://pyco2sys.readthedocs.io/en/latest/validate/)
    varargout{1}=Revelle;
end


function varargout=CalculateAlkParts(pH,pH_scale,Ks,boron_concentration,fluorine_concentration,sulphate_concentration,phosphate_concentration,silicate_concentration,ammonia_concentration,sulphide_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    % ' SUB CalculateAlkParts, version 01.03, 10-10-97, written by Ernie Lewis.
    % ' Inputs: pH, TC, K(), T()
    % ' Outputs: BAlk, OH, PAlk, SiAlk, Hfree, HSO4, HF
    % ' This calculates the various contributions to the alkalinity.
    % ' Though it is coded for H on the total pH scale, for the pH values occuring
    % ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
    % ' negligible) as long as the K Constants are on that scale.
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    
    KWF =KW(selected);
    KP1F=KP1(selected);   KP2F=KP2(selected);   KP3F=KP3(selected);   TPF=phosphate_concentration(selected);
    TSiF=silicate_concentration(selected);   KSiF=KSi(selected);   TNH4F=ammonia_concentration(selected); KNH4F=KNH4(selected);
    TH2SF=sulphide_concentration(selected); KH2SF=KH2S(selected); TBF =boron_concentration(selected);    KBF=KB(selected);
    TSF =sulphate_concentration(selected);    KSF =KS(selected);    TFF =fluorine_concentration(selected);    KFF=KF(selected);
    
    H         = 10.^(-pH(selected));
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH(selected),pH_scale,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree); %' since KF is on the free scale
    
    varargout{1} = BAlk;  varargout{2} = OH; varargout{3} = PAlk;
    varargout{4} = SiAlk; varargout{5} = AmmAlk; varargout{6} = HSAlk;
    varargout{7} = Hfree; varargout{8} = HSO4; varargout{9} = HF;
end


function varargout=CaSolubility(salinity, TempC, TC, pH, Ks, sqrt_salinity, gas_constant, calcium_concentration,which_k1_k2,Pbar,selected)
    global k_perturbation_GLOBAL    % Id of perturbed K
    global Perturb % perturbation
    % '***********************************************************************
    % ' SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
    % ' Inputs: which_k1_k2%, salinity, temperature_in, pressure_in, TCi, pHi, K1, K2
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
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);

    temp_k    = TempC + 273.15;
    log_temp_k = log(temp_k);

    Ca=calcium_concentration(selected);
    Ar=nan(sum(selected),1);
    KCa=nan(sum(selected),1);
    KAr=nan(sum(selected),1);
    TempKx=temp_k;
    logTempKx=log_temp_k;
    sqrSalx=sqrt_salinity;
    Pbarx=Pbar(selected);
    RR = (gas_constant.*temp_k);
    RTx = RR;
    FF=(which_k1_k2(selected)~=6 & which_k1_k2(selected)~=7);
    if any(FF)
    % (below here, selected isn't used, since almost always all rows match the above criterium,
    %  in all other cases the rows will be overwritten later on).
        % CalciteSolubility:
        % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKCa = -171.9065 - 0.077993.*TempKx(FF) + 2839.319./TempKx(FF);
        logKCa = logKCa + 71.595.*logTempKx(FF)./log(10);
        logKCa = logKCa + (-0.77712 + 0.0028426.*TempKx(FF) + 178.34./TempKx(FF)).*sqrSalx(FF);
        logKCa = logKCa - 0.07711.*salinity(FF) + 0.0041249.*sqrSalx(FF).*salinity(FF);
        % '       sd fit = .01 (for salinity part, not part independent of salinity)
        KCa(FF) = 10.^(logKCa);% ' this is in (mol/kg-SW)^2
        % AragoniteSolubility:
        % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKAr = -171.945 - 0.077993.*TempKx(FF) + 2903.293./TempKx(FF);
        logKAr = logKAr + 71.595.*logTempKx(FF)./log(10);
        logKAr = logKAr + (-0.068393 + 0.0017276.*TempKx(FF) + 88.135./TempKx(FF)).*sqrSalx(FF);
        logKAr = logKAr - 0.10018.*salinity(FF) + 0.0059415.*sqrSalx(FF).*salinity(FF);
        % '       sd fit = .009 (for salinity part, not part independent of salinity)
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
    FF=(which_k1_k2(selected)==6 | which_k1_k2(selected)==7);
    if any(FF)
        % *** CalculateKCaforGEOSECS:
        % Ingle et al, Marine Chemistry 1:295-307, 1973 is referenced in
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
        % but the fit is actually from Ingle, Marine Chemistry 3:301-319, 1975)
        KCa(FF) = 0.0000001.*(-34.452 - 39.866.*salinity(FF).^(1./3) +...
            110.21.*log(salinity(FF))./log(10) - 0.0000075752.*TempKx(FF).^2);
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
    % Requested perturbation is passed through global variables k_perturbation_GLOBAL and Perturb
    if (~ isempty(k_perturbation_GLOBAL))
        switch k_perturbation_GLOBAL
            case {'KSPA'}   % solubility Product for Aragonite
                KAr = KAr + Perturb;
            case {'KSPC'}   % for Calcite
                KCa = KCa + Perturb;
            case {'calcium_concentration'}   % for calcium concentration
                Ca  = Ca  + Perturb;
        end
    end
    
    % CalculateOmegasHere:
    H = 10.^(-pH);
    CO3 = TC.*K1(selected).*K2(selected)./(K1(selected).*H + H.*H + K1(selected).*K2(selected));
    varargout{1} = CO3.*Ca./KCa; % OmegaCa, dimensionless
    varargout{2} = CO3.*Ca./KAr; % OmegaAr, dimensionless
end
    
function varargout=FindpHOnAllScales(pH,pH_scale_in,Ks,fluorine_concentration,sulphate_concentration,selected,number_of_points,which_k1_k2,salinity,temp_k)
    % ' SUB FindpHOnAllScales, version 01.02, 01-08-97, written by Ernie Lewis.
    % ' Inputs: pHScale%, pH, K(), T(), fH
    % ' Outputs: pHNBS, pHfree, pHTot, pHSWS
    % ' This takes the pH on the given scale and finds the pH on all scales.
    %  sulphate_concentration = T(3); fluorine_concentration = T(2);
    %  KS = K(6); KF = K(5);% 'these are at the given T, S, P
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = unpack_Ks(Ks);
    fH = calculate_fH(number_of_points,which_k1_k2,salinity,temp_k);

    TSx=sulphate_concentration(selected); KSx=KS(selected); TFx=fluorine_concentration(selected); KFx=KF(selected);fHx=fH(selected);
    FREEtoTOT = (1 + TSx./KSx); % ' pH scale conversion factor
    SWStoTOT  = (1 + TSx./KSx)./(1 + TSx./KSx + TFx./KFx);% ' pH scale conversion factor
    factor=nan(sum(selected),1);
    nF=pH_scale_in(selected)==1;  %'"pHtot"
    factor(nF) = 0;
    nF=pH_scale_in(selected)==2; % '"pHsws"
    factor(nF) = -log(SWStoTOT(nF))./log(0.1);
    nF=pH_scale_in(selected)==3; % '"pHfree"
    factor(nF) = -log(FREEtoTOT(nF))./log(0.1);
    nF=pH_scale_in(selected)==4;  %'"pHNBS"
    factor(nF) = -log(SWStoTOT(nF))./log(0.1) + log(fHx(nF))./log(0.1);
    pHtot  = pH    - factor;    % ' pH comes into this sub on the given scale
    pHNBS  = pHtot - log(SWStoTOT) ./log(0.1) + log(fHx)./log(0.1);
    pHfree = pHtot - log(FREEtoTOT)./log(0.1);
    pHsws  = pHtot - log(SWStoTOT) ./log(0.1);
    varargout{1}=pHtot;
    varargout{2}=pHsws;
    varargout{3}=pHfree;
    varargout{4}=pHNBS;
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
