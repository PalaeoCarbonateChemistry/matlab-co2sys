function VPFac = calculate_VPFac(salinity)
    global temp_k_GLOBAL
    % CalculateVPFac:
    % Weiss, R. selected_GLOBAL., and Price, B. A., Nitrous oxide solubility in water and
    %       seawater, Marine Chemistry 8:347-359, 1980.
    % They fit the data of Goff and Gratch (1946) with the vapor pressure
    %       lowering by sea salt as given by Robinson (1954).
    % This fits the more complicated Goff and Gratch, and Robinson equations
    %       from 273 to 313 deg K and 0 to 40 Sali with a standard error
    %       of .015%, about 5 uatm over this range.
    % This may be on IPTS-29 since they didn't mention the temperature scale,
    %       and the data of Goff and Gratch came before IPTS-48.
    % The references are:
    % Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
    %       to 212 deg selected_GLOBAL, Transactions of the American Society of Heating and
    %       Ventilating Engineers 52:95-122, 1946.
    % Robinson, Journal of the Marine Biological Association of the U. K.
    %       33:449-455, 1954.
    %       This is eq. 10 on p. 350.
    %       This is in atmospheres.
    VPWP = exp(24.4543 - 67.4509.*(100./temp_k_GLOBAL) - 4.8489.*log(temp_k_GLOBAL./100));
    VPCorrWP = exp(-0.000544.*salinity);
    VPSWWP = VPWP.*VPCorrWP;
    VPFac = 1 - VPSWWP; % this assumes 1 atmosphere
end