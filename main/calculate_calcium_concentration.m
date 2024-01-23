function calcium_concentration = calculate_calcium_concentration(salinity,number_of_points,which_k1_k2)
    % Generate empty vectors for holding results
    calcium_concentration = nan(number_of_points,1);    
    
    % CalculateCAL - Total Calcium:
    selected=(which_k1_k2~=6 & which_k1_k2~=7);
        % Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
        % this is .010285.*Sali./35
        calcium_concentration(selected) = 0.02128./40.087.*(salinity(selected)./1.80655);


    selected=(which_k1_k2==6 | which_k1_k2==7);
        % *** CalculateCaforGEOSECS:
        % Culkin, selected, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982) 
        calcium_concentration(selected) = 0.01026.*salinity(selected)./35;
end