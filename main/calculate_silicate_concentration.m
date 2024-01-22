function silicate_concentration = calculate_silicate_concentration(silicate_input,number_of_points)
global which_k1_k2_constants_GLOBAL

    silicate_concentration = nan(number_of_points,1);

    selected=(which_k1_k2_constants_GLOBAL==8 | which_k1_k2_constants_GLOBAL==6);  
    silicate_concentration(selected)  = 0;


    selected=~selected;                         
    silicate_concentration(selected)   = silicate_input(selected)./1e6;
end