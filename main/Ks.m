classdef Ks
    properties
        controls

        k0
        k1
        k2
        kw
        kb
        kf
        ks
        kp1
        kp2
        kp3
        ksi
        kc
        ka
        knh4
        kh2s
    end
    methods
        function self = Ks(temperature_celcius,pressure,composition,pH_scale,which_ks,co2_pressure_correction)
            self.controls = KControls(temperature_celcius,pressure,composition,pH_scale,which_ks,co2_pressure_correction);
        end

        function self = calculate_k0(self)
            self.k0 = KFunctions.calculate_k0(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.co2_pressure_correction,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_k1(self)
            self.k1 = KFunctions.calculate_k1(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_k2(self)
            self.k2 = KFunctions.calculate_k2(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kw(self)
            self.kw = KFunctions.calculate_kw(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kb(self)
            self.kb = KFunctions.calculate_kb(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kf(self)
            self.kf = KFunctions.calculate_kf(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks);
        end
        function self = calculate_ks(self)
            self.ks = KFunctions.calculate_ks(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks);
        end
        function self = calculate_kp1(self)
            self.kp1 = KFunctions.calculate_kp1(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kp2(self)
            self.kp2 = KFunctions.calculate_kp2(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kp3(self)
            self.kp3 = KFunctions.calculate_kp3(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_ksi(self)
            self.ksi = KFunctions.calculate_ksi(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kc(self)
            self.kc = KFunctions.calculate_kc(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_ka(self)
            self.ka = KFunctions.calculate_ka(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_knh4(self)
            self.knh4 = KFunctions.calculate_knh4(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
        function self = calculate_kh2s(self)
            self.kh2s = KFunctions.calculate_kh2s(self.controls.temperature_celcius,self.controls.pressure,self.controls.composition.salinity,self.controls.which_ks,self.controls.pH_scale_conversion);
        end
    
        function self = calculate_all(self)
            self = self.calculate_k0()...
                       .calculate_k1()...
                       .calculate_k2()...
                       .calculate_kw()...
                       .calculate_kb()...
                       .calculate_kf()...
                       .calculate_ks()...
                       .calculate_kp1()...
                       .calculate_kp2()...
                       .calculate_kp3()...
                       .calculate_ksi()...
                       .calculate_kc()...
                       .calculate_ka()...
                       .calculate_knh4()...
                       .calculate_kh2s();
        end

        function [k0,k1,k2,kw,kb,kf,ks,kp1,kp2,kp3,ksi,kc,ka,knh4,kh2s] = unpack_all(self)
            k0 = self.k0;
            k1 = self.k1;
            k2 = self.k2;
            kw = self.kw;
            kb = self.kb;
            kf = self.kf;
            ks = self.ks;
            kp1 = self.kp1;
            kp2 = self.kp2;
            kp3 = self.kp3;
            ksi = self.ksi;
            kc = self.kc;
            ka = self.ka;
            knh4 = self.knh4;
            kh2s = self.kh2s;
        end
        function varargout = unpack_some(self,names)
            for index = 1:numel(names)
                varargout{index} = self.(lower(names(index)));
            end
        end
        function k_map = as_map(self)
            k_map = containers.Map(["K0","K1","K2","KW","KB","KF","KS","KP1","KP2","KP3","KSi","KC","KA","KNH4","KH2S"], ...
                        {self.k0,self.k1,self.k2,self.kw,self.kb,self.kf,self.ks,self.kp1,self.kp2,self.kp3,self.ksi,self.kc,self.ka,self.knh4,self.kh2s});
        end
    
        function ks = select(self,selected)
            if sum(selected)>0
                if numel(selected)==1 || all(selected)
                    ks = self;
                else
                    ks = Ks(self.controls.temperature_celcius(selected),self.controls.pressure(selected),self.controls.composition.select(selected),...
                        self.controls.which_pH_scale(selected),self.controls.which_ks.select(selected),self.controls.co2_pressure_correction(selected));

                    ks.k0 = self.k0(selected);
                    ks.k1 = self.k1(selected);
                    ks.k2 = self.k2(selected);
                    ks.kw = self.kw(selected);
                    ks.kb = self.kb(selected);
                    ks.kf = self.kf(selected);
                    ks.ks = self.ks(selected);
                    ks.kp1 = self.kp1(selected);
                    ks.kp2 = self.kp2(selected);
                    ks.kp3 = self.kp3(selected);
                    ks.kc = self.kc(selected);
                    ks.ka = self.ka(selected);
                    ks.ksi = self.ksi(selected);
                    ks.knh4 = self.knh4(selected);
                    ks.kh2s = self.kh2s(selected);
                end
            else
                ks = NaN;
            end
        end
    end
end