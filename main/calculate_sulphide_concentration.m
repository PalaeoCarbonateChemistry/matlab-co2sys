function sulphide_concentration = calculate_sulphide_concentration(sulphide_input)
global which_k1_k2_constants_GLOBAL number_of_points

    sulphide_concentration = nan(number_of_points,1);

    selected=(which_k1_k2_constants_GLOBAL==8 | which_k1_k2_constants_GLOBAL==6);  
    sulphide_concentration(selected)  = 0;


    selected=~selected;                         
    sulphide_concentration(selected)   = sulphide_input(selected)./1e6;
end