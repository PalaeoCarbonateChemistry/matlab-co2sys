function ammonia_concentration = calculate_ammonia_concentration(ammonia_input)
global which_k1_k2_constants_GLOBAL number_of_points

    ammonia_concentration = nan(number_of_points,1);

    selected=(which_k1_k2_constants_GLOBAL==8 | which_k1_k2_constants_GLOBAL==6);  
    ammonia_concentration(selected)  = 0;


    selected=~selected;                         
    ammonia_concentration(selected)   = ammonia_input(selected)./1e6;
end