%% Defines the class containing tests which use the GLODAP dataset
classdef GlodapTest < matlab.unittest.TestCase
    methods (Static)
        % CO2SYS related methods
        function co2sys = calculate_from_dic_and_alkalinity()
            % calculate_from_dic_and_alkalinity  
            % Calculates carbonate system properties using DIC and alkalinity in GLODAP dataset
            % Returns a matrix describing the carbonate system (as
            % described in the CO2SYS documentation)
            glodap_data = GlodapTest.get_glodap_subset();
            co2sys = CO2SYS(glodap_data.dic,glodap_data.alkalinity,2,1,glodap_data.salinity,glodap_data.temperature,glodap_data.temperature,glodap_data.pressure,glodap_data.pressure,glodap_data.silicate,glodap_data.phosphate,0,0,1,4,1,1,1);
        end
        function co2sys = calculate_from_pH_and_dic()
            % calculate_from_pH_and_dic  
            % Calculates carbonate system properties using pH and DIC in GLODAP dataset
            % Returns a matrix describing the carbonate system (as
            % described in the CO2SYS documentation)
            glodap_data = GlodapTest.get_glodap_subset();
            co2sys = CO2SYS(glodap_data.dic,glodap_data.pH,2,3,glodap_data.salinity,glodap_data.temperature,glodap_data.temperature,glodap_data.pressure,glodap_data.pressure,glodap_data.silicate,glodap_data.phosphate,0,0,1,4,1,1,1);
        end
        function co2sys = calculate_from_pH_and_alkalinity()
            % calculate_from_pH_and_alkalinity 
            % Calculates carbonate system properties using pH and DIC in GLODAP dataset
            % Returns a matrix describing the carbonate system (as
            % described in the CO2SYS documentation)
            glodap_data = GlodapTest.get_glodap_subset();
            co2sys = CO2SYS(glodap_data.pH,glodap_data.alkalinity,3,1,glodap_data.salinity,glodap_data.temperature,glodap_data.temperature,glodap_data.pressure,glodap_data.pressure,glodap_data.silicate,glodap_data.phosphate,0,0,1,4,1,1,1);
        end
        
        % Glodap related methods
        function success = download_glodap_data(inputs)
            % download_glodap_data
            % Downloads a version of glodap v2 data
            % Inputs: string describing version (either
            % "2019","2020","2021","2022")
            % Outputs: flag describing whether acquisition was sucessful
            % Saves file into the /tests/data/ directory
            arguments
                inputs.version = "2022"
            end
            
            try
                links = containers.Map(["2019","2020","2021","2022"], ...
                                     ["https://www.nodc.noaa.gov/archive/arc0133/0186803/1.1/data/0-data/GLODAPv2.2019_Merged_Master_File.mat", ...
                                     "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0210813/GLODAPv2.2020_Merged_Master_File.mat", ...
                                     "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0237935/GLODAPv2.2021_Merged_Master_File.mat", ...
                                     "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0257247/GLODAPv2.2022_Merged_Master_File.mat"]);
                
                websave("./../data/glodap_"+inputs.version,links(inputs.version));
        
                success = true;
            catch
                success = false;
            end
        end
        function data = get_glodap_data(inputs)
            % get_glodap_data
            % Acquires glodap data either by loading from file or
            % downloading if not available
            % Inputs: string describing version (either
            % "2019","2020","2021","2022")
            % Outputs: Variable of loaded data
            arguments
                inputs.version = "2022"
            end
        
            if ~exist("./../data/glodap_"+inputs.version,"file")==2
                downloaded = GlodapTest.download_glodap_data(version=inputs.version);
                if downloaded
                    data = load("./../data/glodap_"+inputs.version);
                else
                    error("Unable to acquire data");
                end
            else
                data = load("./../data/glodap_"+inputs.version);
            end
        end
        function data = get_glodap_subset(inputs)
            % get_glodap_subset
            % Acquires glodap data and translates it into a useful form,
            % with extraneous variables and useless data removed
            % Inputs: string describing version (either
            % "2019","2020","2021","2022")
            % Outputs: A table of valid data for comparison
            arguments
                inputs.version = "2022"
            end
            glodap_data = GlodapTest.get_glodap_data(version=inputs.version);
        
            bottom_depth = glodap_data.G2bottomdepth;
            depth = glodap_data.G2depth;
        
            [aou,aou_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"aou");
            [dic,dic_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"tco2");
            [alkalinity,alkalinity_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"talk");
            [pH,pH_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"phtsinsitutp");
            [pH_standard_conditions,pH_standard_conditions_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"phts25p0");
            
            temperature = glodap_data.G2temperature;
            [salinity,salinity_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"salinity");
            pressure = glodap_data.G2pressure;
        
            [phosphate,phosphate_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"phosphate");
            [silicate,silicate_flag] = GlodapTest.extract_raw_glodap_variable(glodap_data,"silicate");
        
            good_data = aou_flag==2 & dic_flag==2 & alkalinity_flag==2 & pH_flag==2 & pH_standard_conditions_flag==2 & salinity_flag==2 & phosphate_flag==2 & silicate_flag==2;
        
            data = table(aou(good_data), ...
                         dic(good_data), ...
                         alkalinity(good_data), ...
                         pH(good_data), ...
                         pH_standard_conditions(good_data), ...
                         temperature(good_data), ...
                         salinity(good_data), ...
                         pressure(good_data), ...
                         phosphate(good_data), ...
                         silicate(good_data), ...
                         VariableNames=["aou","dic","alkalinity","pH","pH_standard_conditions","temperature","salinity","pressure","phosphate","silicate"]);    
        end
        function generate_glodap_results()
            % generate_glodap_results 
            % Calculates carbonate system properties using three combinations of pH, DIC, and alkalinity in GLODAP dataset
            % Saves the resulting carbonate systems to file
            % To be run when producing datasets for CO2SYS to verify
            % against (i.e. only when a change is made to the code which
            % shoudl change the result).
            co2sys_dic_alkalinity = GlodapTest.calculate_from_dic_and_alkalinity();
            co2sys_pH_alkalinity = GlodapTest.calculate_from_pH_and_alkalinity();
            co2sys_pH_dic = GlodapTest.calculate_from_pH_and_dic();

            save("./../data/glodap_dic_alkalinity.mat","co2sys_dic_alkalinity");
            save("./../data/glodap_pH_alkalinity.mat","co2sys_pH_alkalinity");
            save("./../data/glodap_pH_dic.mat","co2sys_pH_dic");
        end
        function [dic_alkalinity,pH_alkalinity,pH_dic] = load_glodap_data()
            % load_glodap_data 
            % Loads the saved glodap files
            % Returns three matrices for the three combinations of glodap
            % carbonate chemistry parameters (dic + alkalinity, pH +
            % alkalinity, and pH + DIC)
            dic_alkalinity = load("./../data/glodap_dic_alkalinity.mat").co2sys_dic_alkalinity;
            pH_alkalinity = load("./../data/glodap_pH_alkalinity.mat").co2sys_pH_alkalinity;
            pH_dic = load("./../data/glodap_pH_dic.mat").co2sys_pH_dic;
        end
        function [variable,flag] = extract_raw_glodap_variable(glodap_data,name)
            variable = glodap_data.("G2"+name);
            flag = glodap_data.("G2"+name+"f");
        end
    end
    %% Test Method Block
    methods (Test)        
        %% Test data aquisition
        function test_download_data(testCase)
            % test_download_data
            % Tests whether data of the default version can be downloaded
            % sucessfully (warning downloaded ~100MB data).
            downloaded = GlodapTest.download_glodap_data();
            
            testCase.verifyEqual(downloaded,true)
        end
        function compare_glodap_results(testCase)
            % compare_glodap_results
            % Uses pre-existing glodap results to verify identical results
            % with current code version
            [pH_calculated,dic_calculated,alkalinity_calculated] = GlodapTest.load_glodap_data();
            
            co2sys_pH_calculated = GlodapTest.calculate_from_dic_and_alkalinity();
            testCase.verifyEqual(pH_calculated,co2sys_pH_calculated);

            co2sys_dic_calculated = GlodapTest.calculate_from_pH_and_alkalinity();
            testCase.verifyEqual(dic_calculated,co2sys_dic_calculated);

            co2sys_alkalinity_calculated = GlodapTest.calculate_from_pH_and_dic();
            testCase.verifyEqual(alkalinity_calculated,co2sys_alkalinity_calculated);
        end
        
    end
end