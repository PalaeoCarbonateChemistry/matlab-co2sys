%% Defines the class containing tests which check incompatible lengths won't pass validation
classdef LengthTests < matlab.unittest.TestCase
    methods (Static)
    end
    methods (Test)
        function all_scalar(testCase)
            parameter_1 = 2000;
            parameter_2 = 8.1;
            parameter_1_type = 2;
            parameter_2_type = 3;

            salinity_in = 35.0;

            temperature_in = 25.0;
            temperature_out = 25.0;

            pressure_in = 5.0;
            pressure_out= 5.0;

            silicate = 65.5391;
            phosphate = 1.7797;
            ammonia = 0.0;
            sulphide = 0.0;

            pH_scale_in = 1;
            which_k1_k2 = 4;
            which_kso4 = 1;
            which_kf = 1;
            which_boron = 1;

            co2_pressure_correction = 1;

            CO2SYS(parameter_1,parameter_2, ...
                   parameter_1_type,parameter_2_type, ...
                   salinity_in, ...
                   temperature_in,temperature_out, ...
                   pressure_in,pressure_out, ...
                   silicate,phosphate,ammonia,sulphide, ...
                   pH_scale_in, ...
                   which_k1_k2,which_kso4,which_kf, which_boron, ...
                   "co2_press",co2_pressure_correction);
        end
        function none_scalar_mismatch(testCase)
            % none_scalar
            % Puts arrays of mixed lengths (2 and 3 long) into CO2SYS
            % These would previously have passed validation, but should
            % still error. This checks that CO2SYS does emit the expected
            % error.
            parameter_1 = [2000,2100];
            parameter_2 = [8.1,8.2];
            parameter_1_type = [2,2];
            parameter_2_type = [3,3];

            salinity_in = [35.0,35.5,36.0];

            temperature_in = [25.0,26.0,27.0];
            temperature_out = [25.0,26.0,27.0];

            pressure_in = [5.0,5.0];
            pressure_out= [5.0,5.0];

            silicate = [65.5391,65.5391];
            phosphate = [1.7797,1.7797];
            ammonia = [0.0,0.0];
            sulphide = [0.0,0.0];

            pH_scale_in = [1,1];
            which_k1_k2 = [4,4];
            which_kso4 = [1,1];
            which_kf = [1,1];
            which_boron = [1,1];

            co2_pressure_correction = [1,1];

            try
                CO2SYS(parameter_1,parameter_2, ...
                       parameter_1_type,parameter_2_type, ...
                       salinity_in, ...
                       temperature_in,temperature_out, ...
                       pressure_in,pressure_out, ...
                       silicate,phosphate,ammonia,sulphide, ...
                       pH_scale_in, ...
                       which_k1_k2,which_kso4,which_kf, which_boron, ...
                       "co2_press",co2_pressure_correction);
            catch error
                if ~string(error.message).startsWith("Input error - ")
                    error("The mismatching lengths weren't caught by the error checker");
                end
            end
        end
        function none_scalar_match(testCase)
            % none_scalar_match
            % Puts arrays of matched lengths into CO2SYS
            % These should pass validation

            lengths_to_try = [3,10,100];

            for length_to_try = lengths_to_try
                parameter_1 = linspace(2000,2500,length_to_try)';
                parameter_2 = linspace(8.0,8.3,length_to_try)';
                parameter_1_type = 2.*ones(length_to_try,1);
                parameter_2_type = 3.*ones(length_to_try,1);
    
                salinity_in = linspace(30,36,length_to_try)';
    
                temperature_in = linspace(0,40,length_to_try)';
                temperature_out = linspace(10,50,length_to_try)';
    
                pressure_in = linspace(0,500,length_to_try)';
                pressure_out= linspace(0,400,length_to_try)';
    
                silicate = linspace(65,66,length_to_try)';
                phosphate = linspace(1.6,1.9,length_to_try)';
                ammonia = zeros(length_to_try,1);
                sulphide = zeros(length_to_try,1);
    
                pH_scale_in = 1.*ones(length_to_try,1);
                which_k1_k2 = 4.*ones(length_to_try,1);
                which_kso4 = 1.*ones(length_to_try,1);
                which_kf = 1.*ones(length_to_try,1);
                which_boron = 1.*ones(length_to_try,1);
    
                co2_pressure_correction = 1.*ones(length_to_try,1);
    
                try
                    CO2SYS(parameter_1,parameter_2, ...
                           parameter_1_type,parameter_2_type, ...
                           salinity_in, ...
                           temperature_in,temperature_out, ...
                           pressure_in,pressure_out, ...
                           silicate,phosphate,ammonia,sulphide, ...
                           pH_scale_in, ...
                           which_k1_k2,which_kso4,which_kf, which_boron, ...
                           "co2_press",co2_pressure_correction);
                catch error
                    if ~string(error.message).startsWith("Input error - ")
                        error("The mismatching lengths weren't caught by the error checker");
                    end
                end
            end
        end
    end
end