classdef EquilibriumConstants
    properties
        temperature_celcius
        pressure
        salinity
        pH_scale
        which_k1_k2
        which_kf
        which_kso4

        sws_to_tot
    end
    methods
        function self = EquilibriumConstants(temperature_celcius,pressure,salinity,pH_scale,which_k1_k2,which_kso4,which_kf)
            self.temperature_celcius = temperature_celcius;
            self.pressure = pressure;
            self.salinity = salinity;
            self.pH_scale = pH_scale;
            self.which_k1_k2 = which_k1_k2;
            self.which_kf = which_kf;
            self.which_kso4 = which_kso4;
        end
        function k0 = calculate_surface_k0(self)
            k0 = EquilibriumConstantsStatic.calculate_surface_k0(self.temperature_celcius,self.salinity,self.which_k1_k2,self.sws_to_tot);
        end
        function k1 = calculate_surface_k1(self)
            k1 = EquilibriumConstantsStatic.calculate_surface_k1(self.temperature_celcius,self.salinity,self.which_k1_k2,self.sws_to_tot);
        end
        function ks = calculate_surface_ks(self)
            k0 = self.calculate_surface_k0();
            k1 = self.calculate_surface_k1();

            ks = [k0,k1];

            % EquilibriumConstantsStatic.pack_Ks()
        end

    end
end