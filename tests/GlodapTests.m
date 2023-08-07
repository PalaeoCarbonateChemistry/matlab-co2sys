%% Defines the class containing tests which use the GLODAP dataset
classdef GlodapTests < matlab.unittest.TestCase
    methods (Static)
        % CO2SYS related methods
        function co2sys = calculate_from_dic_and_alkalinity(glodap_data)
            % calculate_from_dic_and_alkalinity  
            % Calculates carbonate system properties using DIC and alkalinity in GLODAP dataset
            % Returns a matrix describing the carbonate system (as
            % described in the CO2SYS documentation)
            co2sys = CO2SYS(glodap_data.dic,glodap_data.alkalinity,2,1,glodap_data.salinity,glodap_data.temperature,glodap_data.temperature,glodap_data.pressure,glodap_data.pressure,glodap_data.silicate,glodap_data.phosphate,0,0,1,4,1,1,1);
        end
        function co2sys = calculate_from_pH_and_dic(glodap_data)
            % calculate_from_pH_and_dic  
            % Calculates carbonate system properties using pH and DIC in GLODAP dataset
            % Returns a matrix describing the carbonate system (as
            % described in the CO2SYS documentation)
            co2sys = CO2SYS(glodap_data.dic,glodap_data.pH,2,3,glodap_data.salinity,glodap_data.temperature,glodap_data.temperature,glodap_data.pressure,glodap_data.pressure,glodap_data.silicate,glodap_data.phosphate,0,0,1,4,1,1,1);
        end
        function co2sys = calculate_from_pH_and_alkalinity(glodap_data)
            % calculate_from_pH_and_alkalinity 
            % Calculates carbonate system properties using pH and DIC in GLODAP dataset
            % Returns a matrix describing the carbonate system (as
            % described in the CO2SYS documentation)

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
                
                websave("./big_data/glodap_"+inputs.version,links(inputs.version));
        
                success = true;
            catch
                success = false;
            end
        end
        
        % Get glodap daya
        function data = get_glodap_data_raw(inputs)
            % get_glodap_data_raw
            % Acquires glodap data either by loading from file or
            % downloading if not available
            % Inputs: string describing version (either
            % "2019","2020","2021","2022")
            % Outputs: Variable of loaded data
            arguments
                inputs.version = "2022"
            end
        
            % Check if the file exists
            % If it does fetch that, otherwise download the data
            if ~(exist("./big_data/glodap_"+inputs.version+".mat","file")==2)
                downloaded = GlodapTests.download_glodap_data(version=inputs.version);
                if downloaded
                    data = load("./big_data/glodap_"+inputs.version);
                else
                    error("Unable to acquire data");
                end
            else
                data = load("./big_data/glodap_"+inputs.version);
            end
        end
        function data = get_glodap_data(inputs)
            % get_glodap_data
            % Acquires glodap data and translates it into a useful form,
            % with extraneous variables and useless data removed
            % Inputs: string describing version (either
            % "2019","2020","2021","2022")
            % Outputs: A table of valid data for comparison
            arguments
                inputs.version = "2022"
            end
            % Get raw data
            glodap_data = GlodapTests.get_glodap_data_raw(version=inputs.version);
        
            % Extract relevant parameters (and their flags where relevant)
            bottom_depth = glodap_data.G2bottomdepth;
            depth = glodap_data.G2depth;
        
            [aou,aou_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"aou");
            [dic,dic_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"tco2");
            [alkalinity,alkalinity_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"talk");
            [pH,pH_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"phtsinsitutp");
            [pH_standard_conditions,pH_standard_conditions_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"phts25p0");
            
            temperature = glodap_data.G2temperature;
            [salinity,salinity_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"salinity");
            pressure = glodap_data.G2pressure;
        
            [phosphate,phosphate_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"phosphate");
            [silicate,silicate_flag] = GlodapTests.extract_raw_glodap_variable(glodap_data,"silicate");
        
            % Only use data where the flag for every parameter indicates it
            % is good
            good_data = aou_flag==2 & dic_flag==2 & alkalinity_flag==2 & pH_flag==2 & pH_standard_conditions_flag==2 & salinity_flag==2 & phosphate_flag==2 & silicate_flag==2;
        
            % Compile the data into a table
            data = table(bottom_depth(good_data), ...
                         depth(good_data), ...
                         aou(good_data), ...
                         dic(good_data), ...
                         alkalinity(good_data), ...
                         pH(good_data), ...
                         pH_standard_conditions(good_data), ...
                         temperature(good_data), ...
                         salinity(good_data), ...
                         pressure(good_data), ...
                         phosphate(good_data), ...
                         silicate(good_data), ...
                         VariableNames=["bottom_depth","depth","aou","dic","alkalinity","pH","pH_standard_conditions","temperature","salinity","pressure","phosphate","silicate"]);    
        end
        function glodap_data_subset = get_glodap_data_subset(inputs)
            % get_glodap_data_subset
            % Gets a decimated version of the useful glodap data
            % Inputs: string describing version (either
            % "2019","2020","2021","2022")
            % Outputs: A table of valid data for comparison
            arguments
                inputs.version = "2022";
            end

            % Check if the file exists
            % If it does fetch that, otherwise download the data and create
            % subset
            if ~(exist("./small_data/glodap_subset_"+inputs.version+".mat","file")==2)            
                glodap_data = GlodapTests.get_glodap_data(version=inputs.version);

                % Randomly (but repeatably) decimate the data
                rng(1);
                keep = rand(height(glodap_data),1)<=0.1;
    
                glodap_data_subset = glodap_data(keep,:);
            else
                glodap_data_subset = load("./small_data/glodap_subset_"+inputs.version+".mat").glodap_data_subset;
            end
        end
        
        % Calculate results from glodap inputs
        function generate_glodap_results()
            % generate_glodap_results 
            % Calculates carbonate system properties using three combinations of pH, DIC, and alkalinity in GLODAP dataset
            % Saves the resulting carbonate systems to file
            % To be run when producing datasets for CO2SYS to verify
            % against (i.e. only when a change is made to the code which
            % should change the result).
            glodap_data = GlodapTests.get_glodap_data();

            co2sys_dic_alkalinity = GlodapTests.calculate_from_dic_and_alkalinity(glodap_data);
            co2sys_pH_alkalinity = GlodapTests.calculate_from_pH_and_alkalinity(glodap_data);
            co2sys_pH_dic = GlodapTests.calculate_from_pH_and_dic(glodap_data);

            save("./big_data/glodap_dic_alkalinity.mat","co2sys_dic_alkalinity");
            save("./big_data/glodap_pH_alkalinity.mat","co2sys_pH_alkalinity");
            save("./big_data/glodap_pH_dic.mat","co2sys_pH_dic");
        end
        function generate_glodap_data_subset(inputs)
            arguments
                inputs.version = "2022";
            end
            glodap_data_subset = GlodapTests.get_glodap_data_subset();
            save("./small_data/glodap_subset_"+inputs.version+".mat","glodap_data_subset");
        end
        function generate_glodap_subset_results()
            % generate_glodap_subset_results 
            % Calculates carbonate system properties using three
            % combinations of pH, DIC, and alkalinity in GLODAP subset
            % Saves the resulting carbonate systems to file
            % To be run when producing datasets for CO2SYS to verify
            % against (i.e. only when a change is made to the code which
            % should change the result).
            glodap_data_subset = GlodapTests.get_glodap_data_subset();

            co2sys_dic_alkalinity = GlodapTests.calculate_from_dic_and_alkalinity(glodap_data_subset);
            co2sys_pH_alkalinity = GlodapTests.calculate_from_pH_and_alkalinity(glodap_data_subset);
            co2sys_pH_dic = GlodapTests.calculate_from_pH_and_dic(glodap_data_subset);

            save("./small_data/glodap_subset_dic_alkalinity.mat","co2sys_dic_alkalinity");
            save("./small_data/glodap_subset_pH_alkalinity.mat","co2sys_pH_alkalinity");
            save("./small_data/glodap_subset_pH_dic.mat","co2sys_pH_dic");
        end
        
        % Load results
        function [dic_alkalinity,pH_alkalinity,pH_dic] = load_glodap_results()
            % load_glodap_results 
            % Loads the saved glodap files
            % Returns three matrices for the three combinations of glodap
            % carbonate chemistry parameters (dic + alkalinity, pH +
            % alkalinity, and pH + DIC)
            dic_alkalinity = load("./big_data/glodap_dic_alkalinity.mat").co2sys_dic_alkalinity;
            pH_alkalinity = load("./big_data/glodap_pH_alkalinity.mat").co2sys_pH_alkalinity;
            pH_dic = load("./big_data/glodap_pH_dic.mat").co2sys_pH_dic;
        end
        function [dic_alkalinity,pH_alkalinity,pH_dic] = load_glodap_subset_results()
            % load_glodap_subset_results 
            % Loads the saved decimated glodap files
            % Returns three matrices for the three combinations of glodap
            % carbonate chemistry parameters (dic + alkalinity, pH +
            % alkalinity, and pH + DIC)
            dic_alkalinity = load("./small_data/glodap_subset_dic_alkalinity.mat").co2sys_dic_alkalinity;
            pH_alkalinity = load("./small_data/glodap_subset_pH_alkalinity.mat").co2sys_pH_alkalinity;
            pH_dic = load("./small_data/glodap_subset_pH_dic.mat").co2sys_pH_dic;
        end
        
        % Utility functions
        function [variable,flag] = extract_raw_glodap_variable(glodap_data,name)
            % extract_raw_glodap_variable
            % Extacts a value and its flag from the glodap_data struct
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
            downloaded = GlodapTests.download_glodap_data();
            
            testCase.verifyEqual(downloaded,true)
        end
        function compare_glodap_results(testCase)
            % compare_glodap_results
            % Uses pre-existing glodap results to verify identical results
            % with current code version
            [pH_calculated,dic_calculated,alkalinity_calculated] = GlodapTests.load_glodap_results();

            glodap_data = GlodapTests.get_glodap_data();
            
            co2sys_pH_calculated = GlodapTests.calculate_from_dic_and_alkalinity(glodap_data);
            testCase.verifyEqual(pH_calculated,co2sys_pH_calculated);

            co2sys_dic_calculated = GlodapTests.calculate_from_pH_and_alkalinity(glodap_data);
            testCase.verifyEqual(dic_calculated,co2sys_dic_calculated);

            co2sys_alkalinity_calculated = GlodapTests.calculate_from_pH_and_dic(glodap_data);
            testCase.verifyEqual(alkalinity_calculated,co2sys_alkalinity_calculated);
        end
        function compare_glodap_subset(testCase)
            % compare_glodap_subset
            % Uses pre-existing glodap subset to verify identical results
            % with current code version
            [pH_calculated,dic_calculated,alkalinity_calculated] = GlodapTests.load_glodap_subset_results();

            glodap_data = GlodapTests.get_glodap_data_subset();
            
            co2sys_pH_calculated = GlodapTests.calculate_from_dic_and_alkalinity(glodap_data);
            testCase.verifyEqual(pH_calculated,co2sys_pH_calculated);

            co2sys_dic_calculated = GlodapTests.calculate_from_pH_and_alkalinity(glodap_data);
            testCase.verifyEqual(dic_calculated,co2sys_dic_calculated);

            co2sys_alkalinity_calculated = GlodapTests.calculate_from_pH_and_dic(glodap_data);
            testCase.verifyEqual(alkalinity_calculated,co2sys_alkalinity_calculated);
        end
    end
end