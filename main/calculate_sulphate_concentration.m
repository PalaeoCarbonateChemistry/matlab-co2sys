function sulphate_concentration = calculate_sulphate_concentration()
global salinity_GLOBAL
    % CalculateTS ;
    % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
    % this is .02824.*Sali./35. = .0008067.*Sali
    sulphate_concentration = (0.14./96.062).*(salinity_GLOBAL./1.80655); % in mol/kg-SW    
end