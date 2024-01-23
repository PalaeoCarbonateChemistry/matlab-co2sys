function fH = calculate_fH(number_of_points,which_k1_k2,salinity,temp_k)
    % CalculatefH
    fH = nan(number_of_points,1);
    % Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
    selected=(which_k1_k2==8);
    if any(selected)
        fH(selected) = 1; % this shouldn't occur in the program for this case
    end
    selected=(which_k1_k2==7);
    if any(selected)
        fH(selected) = 1.29 - 0.00204.*  temp_k(selected) + (0.00046 -...
            0.00000148.*temp_k(selected)).*salinity(selected).*salinity(selected);
        % Peng et al, Tellus 39B:439-458, 1987:
        % They reference the GEOSECS report, but round the value
        % given there off so that it is about .008 (1%) lower. It
        % doesn't agree with the check value they give on p. 456.
    end
    selected=(which_k1_k2~=7 & which_k1_k2~=8);
    if any(selected)
        fH(selected) = 1.2948 - 0.002036.*temp_k(selected) + (0.0004607 -...
            0.000001475.*temp_k(selected)).*salinity(selected).^2;
        % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
        % v. 3, 1982 (p. 80);
    end
    
end