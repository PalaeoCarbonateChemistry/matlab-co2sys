function calcium_concentration = calculate_calcium_concentration()
global which_k1_k2_constants_GLOBAL salinity_GLOBAL selected_GLOBAL number_of_points
    % Generate empty vectors for holding results
    calcium_concentration = nan(number_of_points,1);    
    
    % CalculateCAL - Total Calcium:
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7);
        % Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
        % this is .010285.*Sali./35
        calcium_concentration(selected_GLOBAL) = 0.02128./40.087.*(salinity_GLOBAL(selected_GLOBAL)./1.80655);


    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
        % *** CalculateCaforGEOSECS:
        % Culkin, selected_GLOBAL, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982) 
        calcium_concentration(selected_GLOBAL) = 0.01026.*salinity_GLOBAL(selected_GLOBAL)./35;
end