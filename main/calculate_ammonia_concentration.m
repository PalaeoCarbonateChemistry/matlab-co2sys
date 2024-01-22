function ammonia_concentration = calculate_ammonia_concentration(ammonia_input,number_of_points,which_k1_k2)

    ammonia_concentration = nan(number_of_points,1);

    selected=(which_k1_k2==8 | which_k1_k2==6);  
    ammonia_concentration(selected)  = 0;


    selected=~selected;                         
    ammonia_concentration(selected)   = ammonia_input(selected)./1e6;
end