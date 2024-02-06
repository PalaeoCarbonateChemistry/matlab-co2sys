function Ks = calculate_equilibrium_constants(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2)
    KS = EquilibriumConstants.calculate_surface_ks(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2);
    KF = EquilibriumConstants.calculate_surface_kf(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2);
    SWStoTOT  = (1 + sulphate_concentration./KS)./(1 + sulphate_concentration./KS + fluorine_concentration./KF);
    FREEtoTOT =  1 + sulphate_concentration./KS;

    temp_k    = temp_c + 273.15;
    log_temp_k = log(temp_k);

    Ks = EquilibriumConstants.calculate_surface_all(number_of_points,temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,fluorine_concentration,sulphate_concentration,which_kf,which_kso4,which_k1_k2,SWStoTOT,FREEtoTOT);
    [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = EquilibriumConstants.unpack_Ks(Ks);

    fH = calculate_fH(number_of_points,which_k1_k2,salinity,temp_k);

    pressure_correction = EquilibriumConstants.calculate_pressure_correction_all(number_of_points,which_k1_k2,gas_constant,temp_c,pressure_bar,co2_correction);
    [K0_pressure_correction,K1_pressure_correction,K2_pressure_correction,KW_pressure_correction,KB_pressure_correction,KF_pressure_correction,KS_pressure_correction,KP1_pressure_correction,KP2_pressure_correction,KP3_pressure_correction,KSi_pressure_correction,KNH4_pressure_correction,KH2S_pressure_correction] = EquilibriumConstants.unpack_Ks(pressure_correction);
    K0 = K0.*K0_pressure_correction;
    K1 = K1.*K1_pressure_correction;
    K2 = K2.*K2_pressure_correction;
    KW = KW.*KW_pressure_correction;
    KB = KB.*KB_pressure_correction;
    KS = KS.*KS_pressure_correction;
    KF = KF.*KF_pressure_correction;
    KP1 = KP1.*KP1_pressure_correction;
    KP2 = KP2.*KP2_pressure_correction;
    KP3 = KP3.*KP3_pressure_correction;
    KSi = KSi.*KSi_pressure_correction;
    KNH4 = KNH4.*KNH4_pressure_correction;
    KH2S = KH2S.*KH2S_pressure_correction;

    % CorrectpHScaleConversionsForPressure:
    % fH has been assumed to be independent of pressure.
    SWStoTOT  = (1 + sulphate_concentration./KS)./(1 + sulphate_concentration./KS + fluorine_concentration./KF);
    FREEtoTOT =  1 + sulphate_concentration./KS;
    
    %  The values KS and KF are already pressure-corrected, so the pH scale
    %  conversions are now valid at pressure.
    
    % FindpHScaleConversionFactor:
    % this is the scale they will be put on
    pHfactor  = nan(number_of_points,1);
    selected=(pH_scale==1); %Total
    pHfactor(selected) = SWStoTOT(selected);
    selected=(pH_scale==2); %SWS, they are all on this now
    pHfactor(selected) = 1;
    selected=(pH_scale==3); %pHfree
    pHfactor(selected) = SWStoTOT(selected)./FREEtoTOT(selected);
    selected=(pH_scale==4); %pHNBS
    pHfactor(selected) = fH(selected);
    
    % ConvertFromSWSpHScaleToChosenScale:
    K1   = K1.* pHfactor;  K2   = K2.* pHfactor;
    KW   = KW.* pHfactor;  KB   = KB.* pHfactor;
    KP1  = KP1.*pHfactor;  KP2  = KP2.*pHfactor;
    KP3  = KP3.*pHfactor;  KSi  = KSi.*pHfactor;
    KNH4 = KNH4.*pHfactor; KH2S = KH2S.*pHfactor;

    Ks = EquilibriumConstants.pack_Ks(K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S);
end
