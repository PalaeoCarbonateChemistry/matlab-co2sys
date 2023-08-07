%% Defines the class containing tests which use the GLODAP dataset
classdef GridTests < matlab.unittest.TestCase
    methods (Static)
        % Grid generation
        function [temperature_array,salinity_array,pressure_array] = generate_temperature_salinity_pressure_grid()
            % generate_temperature_salinity_pressure_grid
            % Creates a grid of a wide range of temperature, salinity, and
            % pressure, returned as flattened arrays
            
            temperature = -10:2:50;
            salinity = 30:2:40;
            pressure = 0:500:10000;

            [temperature_grid,salinity_grid,pressure_grid] = ndgrid(temperature,salinity,pressure);

            temperature_array = temperature_grid(:);
            salinity_array = salinity_grid(:);
            pressure_array = pressure_grid(:);
        end
        function [dic_array,alkalinity_array] = generate_dic_alkalinity_grid()
            % generate_dic_alkalinity_grid
            % Creates a grid of dic and alkalinity,
            % returned as flattened arrays
            
            dic = 1500:10:2500;
            alkalinity = 1700:10:2700;

            [dic_grid,alkalinity_grid] = ndgrid(dic,alkalinity);

            dic_array = dic_grid(:);
            alkalinity_array = alkalinity_grid(:);
        end
        
        % Calculate methods
        function co2sys = calculate_temperature_salinity_pressure_results()
            % calculate_temperature_salinity_pressure_results
            % Creates grids of temperature, salinity, and pressure, and
            % assumes specific DIC and alkalinity to perform the CO2SYS
            % calculation
            % Returns: a co2sys matrix

            [temperature,salinity,pressure] = GridTests.generate_temperature_salinity_pressure_grid();
            dic = 2000;
            alkalinity = 2300;

            co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature,temperature,pressure,pressure,65.5391,1.7797,0,0,1,4,1,1,1);
        end
        function co2sys = calculate_dic_alkalinity_results()
            % calculate_dic_alkalinity_results
            % Creates grids of dic and alkalinity, and
            % assumes specific temperature, salinity, and pressure, to
            % perform the CO2SYS calculation
            % Returns: a co2sys matrix
            
            [dic,alkalinity] = GridTests.generate_dic_alkalinity_grid();

            temperature = 25;
            salinity = 35;
            pressure = 0;

            co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature,temperature,pressure,pressure,65.5391,1.7797,0,0,1,4,1,1,1);
        end
        
        % Reference creation
        function generate_temperature_salinity_pressure_reference()
            % generate_temperature_salinity_pressure_reference
            % Recalculates the co2sys matrix for the grid defined in grid
            % generation methods, and saves it to
            % temperature_salinity_pressure_grid.mat. This should only be
            % run when CO2SYS is updated in such a way that the results are
            % expected to change.
            co2sys = GridTests.calculate_temperature_salinity_pressure_results();
            save("./small_data/temperature_salinity_pressure_grid.mat","co2sys");
        end
        function generate_dic_alkalinity_reference()
            % generate_dic_alkalinity_reference
            % Recalculates the co2sys matrix for the grid defined in grid
            % generation methods, and saves it to
            % dic_alkalinity_grid.mat. This should only be
            % run when CO2SYS is updated in such a way that the results are
            % expected to change.

            co2sys = GridTests.calculate_dic_alkalinity_results();
            save("./small_data/dic_alkalinity_grid.mat","co2sys");
        end
        
        % Load reference states
        function co2sys_reference = load_temperature_salinity_pressure_reference()
            % load_temperature_salinity_pressure_reference
            % Load the reference state for comparison
            co2sys_reference = load("./small_data/temperature_salinity_pressure_grid.mat").co2sys;
        end
        function co2sys_reference = load_dic_alkalinity_reference()
            % load_dic_alkalinity_reference
            % Load the reference state for comparison
            co2sys_reference = load("./small_data/dic_alkalinity_grid.mat").co2sys;
        end
    end
    %% Test Method Block
    methods (Test)
        function test_temperature_salinity_pressure(testCase)
            % test_temperature_salinity_pressure
            % Loads the reference state for a temperature, pressure, and
            % salinity grid, then uses the current version of CO2SYS to
            % recalculate results on the same grid and ensures parity

            co2sys_reference = GridTests.load_temperature_salinity_pressure_reference();
            co2sys_results = GridTests.calculate_temperature_salinity_pressure_results();

            testCase.verifyEqual(co2sys_reference,co2sys_results);
        end
        function test_dic_alkalinity(testCase)
            % test_dic_alkalinity
            % Loads the reference state for a dic and alkalinity grid, then
            % uses the current version of CO2SYS to
            % recalculate results on the same grid and ensures parity

            co2sys_reference = GridTests.load_dic_alkalinity_reference();
            co2sys_results = GridTests.calculate_dic_alkalinity_results();

            testCase.verifyEqual(co2sys_reference,co2sys_results);
        end
    end
end