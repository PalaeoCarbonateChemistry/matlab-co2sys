%% Defines the class containing tests which are chained
classdef ChainTests < matlab.unittest.TestCase
    methods (Static)

    end
    methods (Test)
        function test_input_combinations(testCase)
            % test_input_combinations
            % This is a round-robin style test that uses the output from
            % each step as input to the next, confirming at each stage that
            % the values produced are close (<0.1%) of the first co2sys
            % calculation.
            % It's particularly useful for coverage - i.e. ensuring all the
            % functions have been used at least once in the tests.
            temperature_in = 10;
            temperature_out = 10;
            salinity = 35;
            pressure_in = 0;
            pressure_out = 0;

            silicate = 65.5391;
            phosphate = 1.7797;
            ammonia = 0;
            sulphide = 0;

            % 1,2 (alkalinity,dic)
            alkalinity = 2300;
            dic = 2000;

            starting_co2sys = CO2SYS(alkalinity,dic,1,2,salinity,temperature_in,temperature_out,pressure_in,pressure_out,silicate,phosphate,ammonia,sulphide,1,4,1,1,1);
            
            % 1,3 (alkalinity,pH)
            alkalinity = starting_co2sys(1);
            pH = starting_co2sys(43);

            co2sys = CO2SYS(alkalinity,pH,1,3,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 1,4 (alkalinity,pCO2)
            alkalinity = co2sys(1);
            pCO2 = co2sys(22);

            co2sys = CO2SYS(alkalinity,pCO2,1,4,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 1,5 (alkalinity,fCO2)
            alkalinity = co2sys(1);
            fCO2 = co2sys(23);

            co2sys = CO2SYS(alkalinity,fCO2,1,5,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 1,6 (alkalinity,HCO3)
            alkalinity = co2sys(1);
            HCO3 = co2sys(24);            
            
            co2sys = CO2SYS(alkalinity,HCO3,1,6,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 1,7 (alkalinity,CO3)
            alkalinity = co2sys(1);
            CO3 = co2sys(25);            
            
            co2sys = CO2SYS(alkalinity,CO3,1,7,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 1,8 (alkalinity,CO2)
            alkalinity = co2sys(1);
            CO2 = co2sys(26);            
            
            co2sys = CO2SYS(alkalinity,CO2,1,8,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 2,3 (dic,pH)
            dic = co2sys(2);
            pH = co2sys(43);            
            
            co2sys = CO2SYS(dic,pH,2,3,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 2,4 (dic,pCO2)
            dic = co2sys(2);
            pCO2 = co2sys(22);            
            
            co2sys = CO2SYS(dic,pCO2,2,4,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 2,5 (dic,fCO2)
            dic = co2sys(2);
            fCO2 = co2sys(23);            
            
            co2sys = CO2SYS(dic,fCO2,2,5,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 2,6 (dic,HCO3)
            dic = co2sys(2);
            HCO3 = co2sys(24);            
            
            co2sys = CO2SYS(dic,HCO3,2,6,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 2,7 (dic,CO3)
            dic = co2sys(2);
            CO3 = co2sys(25);            
            
            co2sys = CO2SYS(dic,CO3,2,7,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 2,8 (dic,CO2)
            dic = co2sys(2);
            CO2 = co2sys(26);
            
            co2sys = CO2SYS(dic,CO2,2,8,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 3,4 (pH,pCO2)
            pH = co2sys(43);
            pCO2 = co2sys(22);
            
            co2sys = CO2SYS(pH,pCO2,3,4,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 3,5 (pH,fCO2)
            pH = co2sys(43);
            fCO2 = co2sys(23);
            
            co2sys = CO2SYS(pH,fCO2,3,5,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 3,6 (pH,HCO3)
            pH = co2sys(43);
            HCO3 = co2sys(24);
            
            co2sys = CO2SYS(pH,HCO3,3,6,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 3,7 (pH,CO3)
            pH = co2sys(43);
            CO3 = co2sys(25);
            
            co2sys = CO2SYS(pH,CO3,3,7,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 3,8 (pH,CO2)
            pH = co2sys(43);
            CO2 = co2sys(26);
            
            co2sys = CO2SYS(pH,CO2,3,8,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 4,6 (pCO2,HCO3)
            pCO2 = co2sys(22);
            HCO3 = co2sys(24);
            
            co2sys = CO2SYS(pCO2,HCO3,4,6,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 4,7 (pCO2,CO3)
            pCO2 = co2sys(22);
            CO3 = co2sys(25);
            
            co2sys = CO2SYS(pCO2,CO3,4,7,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 5,6 (fCO2,HCO3)
            fCO2 = co2sys(23);
            HCO3 = co2sys(24);
            
            co2sys = CO2SYS(fCO2,HCO3,5,6,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 5,7 (fCO2,CO3)
            fCO2 = co2sys(23);
            CO3 = co2sys(25);
            
            co2sys = CO2SYS(fCO2,CO3,5,7,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 6,7 (HCO3,CO3)
            HCO3 = co2sys(24);
            CO3 = co2sys(25);
            
            co2sys = CO2SYS(HCO3,CO3,6,7,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 6,8 (HCO3,CO2)
            HCO3 = co2sys(24);
            CO2 = co2sys(26);
            
            co2sys = CO2SYS(HCO3,CO2,6,8,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);

            % 7,8 (CO3,CO2)
            CO3 = co2sys(25);
            CO2 = co2sys(26);
            
            co2sys = CO2SYS(CO3,CO2,7,8,salinity,temperature_in,temperature_out,pressure_in,pressure_out,65.5391,1.7797,0,0,1,4,1,1,1);
            co2sys(51:52) = starting_co2sys(51:52);
            testCase.verifyEqual(starting_co2sys,co2sys,"RelTol",1e-4);
        end
    end
end