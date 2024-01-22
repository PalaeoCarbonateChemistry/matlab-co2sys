function boron_concentration = calculate_boron_concentration(salinity,number_of_points,which_boron)
    global selected_GLOBAL which_k1_k2_constants_GLOBAL

    boron_concentration = NaN(number_of_points,1);
    % CalculateTB - Total Borate:
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==8); % Pure water case.
    if any(selected_GLOBAL)
        boron_concentration(selected_GLOBAL) = 0;
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL==6 | which_k1_k2_constants_GLOBAL==7);
    if any(selected_GLOBAL)
        boron_concentration(selected_GLOBAL) = 0.0004106.*salinity(selected_GLOBAL)./35; % in mol/kg-SW
        % this is .00001173.*Sali
        % this is about 1% lower than Uppstrom's value
        % Culkin, selected_GLOBAL., in Chemical Oceanography,
        % ed. Riley and Skirrow, 1965:
        % GEOSECS references this, but this value is not explicitly
        % given here
    end
    selected_GLOBAL=(which_k1_k2_constants_GLOBAL~=6 & which_k1_k2_constants_GLOBAL~=7 & which_k1_k2_constants_GLOBAL~=8); % All other cases
    if any(selected_GLOBAL)
	    FF=selected_GLOBAL&(which_boron==1); % If user opted for Uppstrom's values:
	    if any(FF)
	        % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
	        % this is .000416.*Sali./35. = .0000119.*Sali
		    % boron_concentration_GLOBAL(FF) = (0.000232./10.811).*(salinity(FF)./1.80655); % in mol/kg-SW
	        boron_concentration(FF) =  0.0004157.*salinity(FF)./35; % in mol/kg-SW
	    end
	    FF=selected_GLOBAL&(which_boron==2); % If user opted for the Lee et al. values:
	    if any(FF)
		    % Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.	
	 	    % Geochimica Et Cosmochimica Acta 74 (6): 1801-1811.
		    boron_concentration(FF) =  0.0004326.*salinity(FF)./35; % in mol/kg-SW
	    end
    end
end