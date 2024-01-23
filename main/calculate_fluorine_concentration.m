function fluorine_concentration = calculate_fluorine_concentration(salinity)
    % CalculateTF;
    % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
    % this is .000068.*Sali./35. = .00000195.*Sali
    fluorine_concentration = (0.000067./18.998).*(salinity./1.80655); % in mol/kg-SW
end