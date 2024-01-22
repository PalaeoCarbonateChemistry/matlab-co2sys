function fugacity_factor = calculate_fugacity_factor(p_opt,gas_constant,number_of_points,which_k1_k2)
    global temp_k_GLOBAL selected_GLOBAL
    % CalculateFugacityConstants:
    % In previos versions of CO2SYS, the fugacity factor was calculated
    % assuming pressure at one atmosphere, or close to it. Starting with
    % v3.2.1, an option to use in situ pressure is provided.
    %       Weiss, R. selected_GLOBAL., Marine Chemistry 2:203-215, 1974.
    %       Delta and B in cm3/mol
    fugacity_factor = ones(number_of_points,1);
    Delta = (57.7 - 0.118.*temp_k_GLOBAL);
    b = -1636.75 + 12.0408.*temp_k_GLOBAL - 0.0327957.*temp_k_GLOBAL.^2 + 3.16528.*0.00001.*temp_k_GLOBAL.^3;
    % For a mixture of CO2 and air at in situ pressure;
    xc2 = 1; % assumed to be 1, though not strictly correct (xc2 = [1-xCO2]^2)
    P1atm = 1.01325; % atmospheric pressure in bar
    if p_opt == 0
        fugacity_factor = exp((b + 2.*xc2.*Delta).*P1atm./(gas_constant.*temp_k_GLOBAL)); % FugFac at 1 atm
    elseif p_opt == 1
        fugacity_factor = exp((b + 2.*xc2.*Delta).*(P1atm+Pbar)./(gas_constant.*temp_k_GLOBAL)); % FugFac at in situ pressure
    else
        disp('co2_press must be set to either 0 or 1'); % Display error message
    end
    selected_GLOBAL=(which_k1_k2==6 | which_k1_k2==7); % GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
    fugacity_factor(selected_GLOBAL) = 1;
end