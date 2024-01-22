function phosphate_concentration = calculate_phosphate_concentration(phosphate_input)
global which_k1_k2_constants_GLOBAL number_of_points

    phosphate_concentration = nan(number_of_points,1);

    selected=(which_k1_k2_constants_GLOBAL==8 | which_k1_k2_constants_GLOBAL==6);  
    phosphate_concentration(selected)  = 0;


    selected=~selected;                         
    phosphate_concentration(selected)   = phosphate_input(selected)./1e6;
end