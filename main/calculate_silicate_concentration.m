function silicate_concentration = calculate_silicate_concentration(silicate_input,number_of_points,which_k1_k2)
    silicate_concentration = nan(number_of_points,1);

    selected=(which_k1_k2==8 | which_k1_k2==6);  
    silicate_concentration(selected)  = 0;


    selected=~selected;                         
    silicate_concentration(selected)   = silicate_input(selected)./1e6;
end