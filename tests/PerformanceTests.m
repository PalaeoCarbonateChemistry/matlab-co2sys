%% Defines the class containing tests which use the GLODAP dataset
classdef PerformanceTests < matlab.unittest.TestCase
    methods (Static)

    end
    methods (Test)
        function test_array_temperature_salinity_pressure(testCase)
            number_of_points = 100000;
            [temperature,salinity,pressure] = RandomTests.generate_temperature_salinity_pressure_arrays(number_of_points);

            dic = 2000;
            alkalinity = 2300;

            number_of_iterations = 100;
            time_taken = NaN(number_of_iterations,1);
    
            for iteration = 1:number_of_iterations
                tic
                co2sys = CO2SYS(dic,alkalinity,2,1,salinity,temperature,temperature,pressure,pressure,65.5391,1.7797,0,0,1,4,1,1,1);
                time_taken(iteration) = toc;
            end

            average_time_taken = median(time_taken);
        end
        function test_iterative_temperature_salinity_pressure(testCase)
            number_of_points = 100;
            [temperature,salinity,pressure] = RandomTests.generate_temperature_salinity_pressure_arrays(number_of_points);

            dic = 2000;
            alkalinity = 2300;

            number_of_iterations = 100;
            time_taken = NaN(number_of_iterations,1);
    
            for iteration = 1:number_of_iterations
                tic
                for point_number = 1:number_of_points
                    local_temperature = temperature(iteration);
                    local_salinity = salinity(iteration);
                    local_pressure = pressure(iteration);

                    co2sys = CO2SYS(dic,alkalinity,2,1,local_salinity,local_temperature,local_temperature,local_pressure,local_pressure,65.5391,1.7797,0,0,1,4,1,1,1);
                end
                time_taken(iteration) = toc;
            end

            average_time_taken = median(time_taken);
        end
    end
end