%% Defines the class containing tests for coverage of EquilibriumConstantsStatic and EquilibriumConstants
classdef EquilibriumConstantsTests < matlab.unittest.TestCase
     methods (Static)
         function calculated = calculate_all_static(inputs)
             temp_c = inputs.temperature_in;
             pressure_bar = inputs.pressure_in;
                
             salinity = inputs.salinity;

             pH_scale = inputs.pH_scale_in;

             silicate = inputs.silicate;
             phosphate = inputs.phosphate;
             ammonia = inputs.ammonia;
             sulphide = inputs.sulphide;
                
             which_k1_k2 = inputs.which_k1k2;
             which_kso4 = inputs.which_kso4;
             which_kf = inputs.which_kf;
             which_ks = WhichKs(which_k1_k2,which_kso4,which_kf);
             which_boron = inputs.which_boron;
             co2_correction = inputs.co2_pressure_correction;

             gas_constant = 83.14462618; % ml bar-1 K-1 mol-1,


             composition = Composition(salinity)...
                            .set_silicate_concentration(silicate/1e6)...
                            .set_phosphate_concentration(phosphate/1e6)...
                            .set_ammonia_concentration(ammonia/1e6)...
                            .set_sulphide_concentration(sulphide/1e6)...
                            .estimate_all_from_salinity(which_boron)...
                            .remove_freshwater_species(which_ks)...
                            .adjust_geosecs_species(which_ks)...
                            .calculate_peng_correction(which_ks);


             calculated = EquilibriumConstantsStatic.calculate_all(temp_c,pressure_bar,salinity,pH_scale,co2_correction,gas_constant,composition,which_ks);
         end
         function calculated = calculate_individual_static(inputs)
             temp_c = inputs.temperature_in;
             pressure_bar = inputs.pressure_in;
                
             salinity = inputs.salinity;
                
             which_k1_k2 = inputs.which_k1k2;
             which_kso4 = inputs.which_kso4;
             which_kf = inputs.which_kf;
             which_ks = WhichKs(which_k1_k2,which_kso4,which_kf);
             which_boron = inputs.which_boron;
             co2_correction = inputs.co2_pressure_correction;

             pH_scale = inputs.pH_scale_in;

             silicate = inputs.silicate;
             phosphate = inputs.phosphate;
             ammonia = inputs.ammonia;
             sulphide = inputs.sulphide;

             gas_constant = 83.14462618; % ml bar-1 K-1 mol-1,
             
             composition = Composition(salinity)...
                            .set_silicate_concentration(silicate/1e6)...
                            .set_phosphate_concentration(phosphate/1e6)...
                            .set_ammonia_concentration(ammonia/1e6)...
                            .set_sulphide_concentration(sulphide/1e6)...
                            .estimate_all_from_salinity(which_boron)...
                            .remove_freshwater_species(which_ks)...
                            .adjust_geosecs_species(which_ks)...
                            .calculate_peng_correction(which_ks);


             ks = EquilibriumConstantsStatic.calculate_surface_ks(temp_c,salinity,which_ks);
             kf = EquilibriumConstantsStatic.calculate_surface_kf(temp_c,salinity,which_ks);
            
             seawater_to_total = EquilibriumConstantsStatic.calculate_seawater_to_total(composition,ks,kf);
             free_to_total = EquilibriumConstantsStatic.calculate_free_to_total(composition,ks);

             ks_deep = EquilibriumConstantsStatic.calculate_ks(temp_c,pressure_bar,salinity,gas_constant,which_ks);
             kf_deep = EquilibriumConstantsStatic.calculate_kf(temp_c,pressure_bar,salinity,gas_constant,which_ks);
            
             seawater_to_total_deep = EquilibriumConstantsStatic.calculate_seawater_to_total(composition,ks_deep,kf_deep);
             free_to_total_deep = EquilibriumConstantsStatic.calculate_free_to_total(composition,ks_deep);

             k0_individual = EquilibriumConstantsStatic.calculate_k0(temp_c,pressure_bar,salinity,co2_correction,gas_constant,which_ks,seawater_to_total);
             k1_individual = EquilibriumConstantsStatic.calculate_k1(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             k2_individual = EquilibriumConstantsStatic.calculate_k2(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             kw_individual = EquilibriumConstantsStatic.calculate_kw(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             kb_individual = EquilibriumConstantsStatic.calculate_kb(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             kp1_individual = EquilibriumConstantsStatic.calculate_kp1(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             kp2_individual = EquilibriumConstantsStatic.calculate_kp2(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             kp3_individual = EquilibriumConstantsStatic.calculate_kp3(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             ksi_individual = EquilibriumConstantsStatic.calculate_ksi(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             knh4_individual = EquilibriumConstantsStatic.calculate_knh4(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);
             kh2s_individual = EquilibriumConstantsStatic.calculate_kh2s(temp_c,pressure_bar,salinity,gas_constant,which_ks,seawater_to_total,seawater_to_total_deep,free_to_total,free_to_total_deep,pH_scale);

             calculated = EquilibriumConstantsStatic.pack_Ks(k0_individual,k1_individual,k2_individual,kw_individual,kb_individual,kf_deep,ks_deep,kp1_individual,kp2_individual,kp3_individual,ksi_individual,knh4_individual,kh2s_individual);
         end
         function calculated = calculate_all_bound(inputs)
             temp_c = inputs.temperature_in;
             pressure_bar = inputs.pressure_in;
                
             salinity = inputs.salinity;

             pH_scale = inputs.pH_scale_in;
                
             which_k1_k2 = inputs.which_k1k2;
             which_kso4 = inputs.which_kso4;
             which_kf = inputs.which_kf;
             co2_correction = inputs.co2_pressure_correction;

             gas_constant = 83.14462618; % ml bar-1 K-1 mol-1,
             sulphate_concentration = calculate_sulphate_concentration(salinity);
             fluorine_concentration = calculate_fluorine_concentration(salinity);

             equilibrium_constants = EquilibriumConstants(temp_c,pressure_bar,salinity,pH_scale,which_k1_k2,which_kso4,which_kf);
             calculated = equilibrium_constants.calculate();
         end
     end
     methods (Test)
         function test_calculate_all_vs_calculate_individual(testCase)
            [co2sys_inputs,~] = RandomTests.load_temperature_salinity_pressure_reference();
            all = EquilibriumConstantsTests.calculate_all_static(co2sys_inputs);
            individual = EquilibriumConstantsTests.calculate_individual_static(co2sys_inputs);

            for key = all.keys()
                testCase.verifyEqual(all(key{1}),individual(key{1}),"RelTol",1e-9);
            end
        end
     end
end