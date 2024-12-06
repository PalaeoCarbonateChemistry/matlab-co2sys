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

             composition = Composition(salinity)...
                            .set_silicate_concentration(silicate/1e6)...
                            .set_phosphate_concentration(phosphate/1e6)...
                            .set_ammonia_concentration(ammonia/1e6)...
                            .set_sulphide_concentration(sulphide/1e6)...
                            .estimate_all_from_salinity(which_boron)...
                            .remove_freshwater_species(which_ks)...
                            .adjust_geosecs_species(which_ks)...
                            .calculate_peng_correction(which_ks);


             calculated = KFunctions.calculate_all(temp_c,pressure_bar,pH_scale,co2_correction,composition,which_ks);
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

             which_pH_scale = inputs.pH_scale_in;

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

            pH_scale_conversion = [pHScale(which_pH_scale,composition,temp_c,composition.salinity,which_ks,0.0),...
                                   pHScale(which_pH_scale,composition,temp_c,composition.salinity,which_ks,pressure_bar)];
            
             % seawater_to_total = KFunctions.calculate_seawater_to_total(composition,ks,kf);
             % free_to_total = KFunctions.calculate_free_to_total(composition,ks);

             % ks_deep = KFunctions.calculate_ks(temp_c,pressure_bar,salinity,gas_constant,which_ks);
             % kf_deep = KFunctions.calculate_kf(temp_c,pressure_bar,salinity,gas_constant,which_ks);
             % 
             % seawater_to_total_deep = KFunctions.calculate_seawater_to_total(composition,ks_deep,kf_deep);
             % free_to_total_deep = KFunctions.calculate_free_to_total(composition,ks_deep);

             k0_individual = KFunctions.calculate_k0(temp_c,pressure_bar,salinity,co2_correction,which_ks,pH_scale_conversion);
             k1_individual = KFunctions.calculate_k1(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             k2_individual = KFunctions.calculate_k2(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             kw_individual = KFunctions.calculate_kw(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             kb_individual = KFunctions.calculate_kb(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             kp1_individual = KFunctions.calculate_kp1(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             kp2_individual = KFunctions.calculate_kp2(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             kp3_individual = KFunctions.calculate_kp3(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             ksi_individual = KFunctions.calculate_ksi(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             
             kc_individual = KFunctions.calculate_kc(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             ka_individual = KFunctions.calculate_ka(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             
             knh4_individual = KFunctions.calculate_knh4(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);
             kh2s_individual = KFunctions.calculate_kh2s(temp_c,pressure_bar,salinity,which_ks,pH_scale_conversion);

             ks_individual = KFunctions.calculate_ks(temp_c,pressure_bar,salinity,which_ks);
             kf_individual = KFunctions.calculate_kf(temp_c,pressure_bar,salinity,which_ks);

             calculated = KFunctions.pack_Ks(k0_individual,k1_individual,k2_individual,kw_individual,kb_individual,kf_individual,ks_individual,kp1_individual,kp2_individual,kp3_individual,ksi_individual,kc_individual,ka_individual,knh4_individual,kh2s_individual);
         end
         function calculated = calculate_all_bound(inputs)
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
             co2_pressure_correction = inputs.co2_pressure_correction;

             composition = Composition(salinity)...
                            .set_silicate_concentration(silicate/1e6)...
                            .set_phosphate_concentration(phosphate/1e6)...
                            .set_ammonia_concentration(ammonia/1e6)...
                            .set_sulphide_concentration(sulphide/1e6)...
                            .estimate_all_from_salinity(which_boron)...
                            .remove_freshwater_species(which_ks)...
                            .adjust_geosecs_species(which_ks)...
                            .calculate_peng_correction(which_ks);

             calculated = Ks(temp_c,pressure_bar,composition,pH_scale,which_ks,co2_pressure_correction)...
                              .calculate_all()...
                              .as_map();
         end
     end
     methods (Test)
         function test_calculate_all_vs_calculate_individual_static(testCase)
            [co2sys_inputs,~] = RandomTests.load_temperature_salinity_pressure_reference();
            all = EquilibriumConstantsTests.calculate_all_static(co2sys_inputs);
            individual = EquilibriumConstantsTests.calculate_individual_static(co2sys_inputs);

            for key = all.keys()
                testCase.verifyEqual(all(key{1}),individual(key{1}),"RelTol",1e-9);
            end
         end
         function test_static_vs_bound(testCase)
            [co2sys_inputs,~] = RandomTests.load_temperature_salinity_pressure_reference();
            static = EquilibriumConstantsTests.calculate_all_static(co2sys_inputs);
            bound = EquilibriumConstantsTests.calculate_all_bound(co2sys_inputs);

            for key = static.keys()
                testCase.verifyEqual(static(key{1}),bound(key{1}),"RelTol",1e-9);
            end
         end
     end
end