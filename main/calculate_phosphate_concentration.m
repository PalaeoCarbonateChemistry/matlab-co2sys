function phosphate_concentration = calculate_phosphate_concentration(phosphate_input,number_of_points,which_k1_k2)
    phosphate_concentration = nan(number_of_points,1);

    selected=(which_k1_k2==8 | which_k1_k2==6);  
    phosphate_concentration(selected)  = 0;


    selected=~selected;                         
    phosphate_concentration(selected)   = phosphate_input(selected)./1e6;
end