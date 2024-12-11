%% Defines the class containing tests which use the GLODAP dataset
classdef RandomTests < matlab.unittest.TestCase
    methods (Static)
        % Utility functions
        function scaled_random = random_between(lower,upper,number_of_points)
            % random_between
            % Returns a number of points between the lower and upper bound
            range = upper-lower;
            random = rand(number_of_points,1);
            scaled_random = (random*range)+lower;
        end
        
        % Array generation
        function [temperature,salinity,pressure] = generate_temperature_salinity_pressure_arrays(number_of_points)
            % generate_temperature_salinity_pressure_arrays
            % Creates arrays for temperature, salinity, and pressure, by
            % random uniform sampling between wide endpoints
            temperature = RandomTests.random_between(-10,50,number_of_points);
            salinity = RandomTests.random_between(0,50,number_of_points);
            pressure = RandomTests.random_between(0,10000,number_of_points);
        end
        function [dic,alkalinity] = generate_dic_alkalinity_arrays(number_of_points)
            % generate_dic_alkalinity_arrays
            % Creates arrays for DIC and alkalinity by
            % random uniform sampling between wide endpoints
            dic = RandomTests.random_between(0,10000,number_of_points);
            alkalinity = RandomTests.random_between(0,10000,number_of_points);
        end        
        
        % Reference creation
        function generate_temperature_salinity_pressure_reference()
            % generate_temperature_salinity_pressure_reference
            % Creates the reference state for 10,000 random combinations of
            % temperature, salinity, and pressure - as well as randomised
            % other inputs (silicate, phosphate, ammonia, sulphide,
            % which_k1k2, which_kso4). DIC and alkalinity are set to
            % approximate modern values and are constant throughout.
            % Input values are saved to a table for easy loading and
            % reprocessing.
            % In theory the generation of random numbers should be
            % repeatable because the seed is set...but this is not certain
            % Saves results to temperature_salinity_pressure_random.mat
            number_of_points = 10000;
            rng(1);

            [temperature_in,salinity,pressure_in] = RandomTests.generate_temperature_salinity_pressure_arrays(number_of_points);
            [temperature_out,~,pressure_out] = RandomTests.generate_temperature_salinity_pressure_arrays(number_of_points);
            
            silicate = RandomTests.random_between(60,70,number_of_points);
            phosphate = RandomTests.random_between(1,5,number_of_points);
            ammonia = round(RandomTests.random_between(0,1,number_of_points),0);
            sulphide = round(RandomTests.random_between(0,1,number_of_points),0);

            pH_scale_in = round(RandomTests.random_between(0,4.5,number_of_points),0);
            
            which_k1k2 = floor(RandomTests.random_between(1,17.5,number_of_points));
            which_kso4 = floor(RandomTests.random_between(1,3.5,number_of_points));
            which_kf = floor(RandomTests.random_between(1,2.5,number_of_points));
            which_boron = floor(RandomTests.random_between(1,2.5,number_of_points));
            co2_pressure_correction = floor(RandomTests.random_between(0,1.5,number_of_points));
            
            dic = repelem(2000,number_of_points,1);
            alkalinity = repelem(2300,number_of_points,1);

            inputs = table(temperature_in,temperature_out, ...
                           salinity, ...
                           pressure_in,pressure_out, ...
                           silicate,phosphate,ammonia,sulphide, ...
                           pH_scale_in,co2_pressure_correction, ...
                           which_k1k2,which_kso4,which_kf,which_boron, ...
                           dic,alkalinity);

            co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature_in,temperature_out,pressure_in,pressure_out,silicate,phosphate,ammonia,sulphide,pH_scale_in,which_k1k2,which_kso4,which_kf,which_boron,co2_pressure_correction);

            writetable(inputs,"./small_data/temperature_salinity_pressure_random_inputs.csv");
            save("./small_data/temperature_salinity_pressure_random.mat","co2sys");
        end
        function generate_dic_alkalinity_reference()
            % generate_dic_alkalinity_reference
            % Creates the reference state for 10,000 random combinations of
            % DIC and alkalinity - as well as randomised
            % other inputs (silicate, phosphate, ammonia, sulphide,
            % which_k1k2, which_kso4). Temperature, salinity, and pressure are set to
            % approximate modern values and are constant throughout.
            % Input values are saved to a table for easy loading and
            % reprocessing.
            % In theory the generation of random numbers should be
            % repeatable because the seed is set...but this is not certain
            % Saves results to dic_alkalinity_random.mat

            number_of_points = 10000;
            rng(1);

            [dic,alkalinity] = RandomTests.generate_dic_alkalinity_arrays(number_of_points);
            
            temperature_in = repelem(25,number_of_points,1);
            temperature_out = repelem(25,number_of_points,1);
            salinity = repelem(35,number_of_points,1);
            pressure_in = repelem(0,number_of_points,1);
            pressure_out = repelem(0,number_of_points,1);

            silicate = RandomTests.random_between(60,70,number_of_points);
            phosphate = RandomTests.random_between(1,5,number_of_points);
            ammonia = round(RandomTests.random_between(0,1,number_of_points),0);
            sulphide = round(RandomTests.random_between(0,1,number_of_points),0);
            
            which_k1k2 = floor(RandomTests.random_between(1,17.5,number_of_points));
            which_kso4 = floor(RandomTests.random_between(1,3.5,number_of_points));

            inputs = table(temperature_in,temperature_out, ...
                           salinity, ...
                           pressure_in,pressure_out, ...
                           silicate,phosphate,ammonia,sulphide, ...
                           which_k1k2,which_kso4, ...
                           dic,alkalinity);

            co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature_in,temperature_out,pressure_in,pressure_out,silicate,phosphate,ammonia,sulphide,1,which_k1k2,which_kso4,1,1);

            writetable(inputs,"./small_data/dic_alkalinity_random_inputs.csv");
            save("./small_data/dic_alkalinity_random.mat","co2sys");
        end
        
        % Reference loading
        function [co2sys_inputs,co2sys_reference] = load_temperature_salinity_pressure_reference()
            % load_temperature_salinity_pressure_reference
            % Loads the reference condition for the random temperature,
            % salinity, and pressure combinations.
            % Returns: Table of inputs to co2sys
            %          The reference co2sys matrix
            co2sys_inputs = readtable("./small_data/temperature_salinity_pressure_random_inputs.csv");
            co2sys_reference = load("./small_data/temperature_salinity_pressure_random.mat").co2sys;
        end
        function [co2sys_inputs,co2sys_reference] = load_dic_alkalinity_reference()
            % load_dic_alkalinity_reference
            % Loads the reference condition for the random DIC and
            % alkalinity combinations.
            % Returns: Table of inputs to co2sys
            %          The reference co2sys matrix
            co2sys_inputs = readtable("./small_data/dic_alkalinity_random_inputs.csv");
            co2sys_reference = load("./small_data/dic_alkalinity_random.mat").co2sys;
        end        
        
        % Result calculation
        function co2sys = calculate_temperature_salinity_pressure_results(inputs)
            % calculate_temperature_pressure_salinity_results
            % Given a table of inputs (loaded from file), unpacks the
            % arguments and runs CO2SYS.
            temperature_in = inputs.temperature_in;
            temperature_out = inputs.temperature_out;

            salinity = inputs.salinity;

            pressure_in = inputs.pressure_in;
            pressure_out = inputs.pressure_out;

            silicate = inputs.silicate;
            phosphate = inputs.phosphate;
            ammonia = inputs.ammonia;
            sulphide = inputs.sulphide;

            pH_scale_in = inputs.pH_scale_in;

            which_k1k2 = inputs.which_k1k2;
            which_kso4 = inputs.which_kso4;
            which_kf = inputs.which_kf;
            which_boron = inputs.which_boron;
            co2_pressure_correction = inputs.co2_pressure_correction;

            dic = inputs.dic;
            alkalinity = inputs.alkalinity;

            co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature_in,temperature_out,pressure_in,pressure_out,silicate,phosphate,ammonia,sulphide,pH_scale_in,which_k1k2,which_kso4,which_kf,which_boron,co2_pressure_correction);
        end
        function co2sys = calculate_dic_alkalinity_results(inputs)
            % calculate_dic_alkalinity_results
            % Given a table of inputs (loaded from file), unpacks the
            % arguments and runs CO2SYS.

            temperature_in = inputs.temperature_in;
            temperature_out = inputs.temperature_out;

            salinity = inputs.salinity;

            pressure_in = inputs.pressure_in;
            pressure_out = inputs.pressure_out;

            silicate = inputs.silicate;
            phosphate = inputs.phosphate;
            ammonia = inputs.ammonia;
            sulphide = inputs.sulphide;

            which_k1k2 = inputs.which_k1k2;
            which_kso4 = inputs.which_kso4;

            dic = inputs.dic;
            alkalinity = inputs.alkalinity;

            co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature_in,temperature_out,pressure_in,pressure_out,silicate,phosphate,ammonia,sulphide,1,which_k1k2,which_kso4,1,1);
        end
    end
    %% Test Method Block
    methods (Test)
        function test_temperature_salinity_pressure(testCase)
            % test_temperature_salinity_pressure
            % Loads the inputs and reference condition from file then
            % calculates the results of CO2SYS using the current version
            % with the same inputs.
            % Verifies these are equal at the relative level of 1e-9

            [co2sys_inputs,co2sys_reference] = RandomTests.load_temperature_salinity_pressure_reference();
            co2sys_results = RandomTests.calculate_temperature_salinity_pressure_results(co2sys_inputs);

            testCase.verifyEqual(co2sys_reference,co2sys_results,"RelTol",1e-9);
        end
        function test_dic_alkalinity(testCase)
            % test_dic_alkalinity
            % Loads the inputs and reference condition from file then
            % calculates the results of CO2SYS using the current version
            % with the same inputs.
            % Verifies these are equal at the relative level of 1e-9

            [co2sys_inputs,co2sys_reference] = RandomTests.load_dic_alkalinity_reference();
            co2sys_results = RandomTests.calculate_dic_alkalinity_results(co2sys_inputs);

            testCase.verifyEqual(co2sys_reference,co2sys_results,"RelTol",1e-9);
        end
    end
end