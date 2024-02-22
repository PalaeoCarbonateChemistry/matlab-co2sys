classdef EquilibriumConstants
    properties
        temperature_celcius
        pressure
        composition

        which_pH_scale
        which_ks

        pH_scale_conversion
        co2_pressure_correction

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
        knh4
        kh2s
    end
    methods
        function self = EquilibriumConstants(temperature_celcius,pressure,composition,pH_scale,which_ks,co2_pressure_correction)
            self.temperature_celcius = temperature_celcius;
            self.pressure = pressure;
            self.composition = composition;
            self.which_pH_scale = pH_scale;
            self.which_ks = which_ks;
            self.co2_pressure_correction = co2_pressure_correction;

            self.pH_scale_conversion = [pHScale(pH_scale,composition,temperature_celcius,composition.salinity,which_ks,0.0),...
                                        pHScale(pH_scale,composition,temperature_celcius,composition.salinity,which_ks,pressure)];
        end
        % function surface_k0 = calculate_surface_k0(self)
        %     surface_k0 = EquilibriumConstantsStatic.calculate_surface_k0(self.temperature_celcius,self.salinity,self.which_ks,self.seawater_to_total);
        % end
        % function self = calculate_surface_ks(self)
        %     self = self.calculate_surface_k0()...;
        %                .calculate_surface_k1();
        % end

        function self = calculate_k0(self)
            self.k0 = EquilibriumConstantsStatic.calculate_k0(self.temperature_celcius,self.pressure,self.composition.salinity,self.co2_pressure_correction,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_k1(self)
            self.k1 = EquilibriumConstantsStatic.calculate_k1(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_k2(self)
            self.k2 = EquilibriumConstantsStatic.calculate_k2(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_kw(self)
            self.kw = EquilibriumConstantsStatic.calculate_kw(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_kb(self)
            self.kb = EquilibriumConstantsStatic.calculate_kb(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_kf(self)
            self.kf = EquilibriumConstantsStatic.calculate_kf(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks);
        end
        function self = calculate_ks(self)
            self.ks = EquilibriumConstantsStatic.calculate_ks(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks);
        end
        function self = calculate_kp1(self)
            self.kp1 = EquilibriumConstantsStatic.calculate_kp1(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_kp2(self)
            self.kp2 = EquilibriumConstantsStatic.calculate_kp2(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_kp3(self)
            self.kp3 = EquilibriumConstantsStatic.calculate_kp3(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_ksi(self)
            self.ksi = EquilibriumConstantsStatic.calculate_ksi(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_knh4(self)
            self.knh4 = EquilibriumConstantsStatic.calculate_knh4(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
        end
        function self = calculate_kh2s(self)
            self.kh2s = EquilibriumConstantsStatic.calculate_kh2s(self.temperature_celcius,self.pressure,self.composition.salinity,self.which_ks,self.pH_scale_conversion);
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
                       .calculate_knh4()...
                       .calculate_kh2s();
        end

        function ks = unpack_ks(self)
            ks = [self.k0,...
                  self.k1,...
                  self.k2,...
                  self.kw,...
                  self.kb,...
                  self.kf,...
                  self.ks,...
                  self.kp1,...
                  self.kp2,...
                  self.kp3,...
                  self.ksi,...
                  self.knh4,...
                  self.kh2s];
        end
        function k_map = as_map(self)
            k_map = containers.Map(["K0","K1","K2","KW","KB","KF","KS","KP1","KP2","KP3","KSi","KNH4","KH2S"], ...
                        {self.k0,self.k1,self.k2,self.kw,self.kb,self.kf,self.ks,self.kp1,self.kp2,self.kp3,self.ksi,self.knh4,self.kh2s});
        end
    end
end