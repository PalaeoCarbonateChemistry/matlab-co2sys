classdef KControls
    properties
        temperature_celcius
        pressure
        composition

        which_pH_scale
        which_ks

        pH_scale_conversion
        co2_pressure_correction
    end
    methods
        function self = KControls(temperature_celcius,pressure,composition,pH_scale,which_ks,co2_pressure_correction)
            self.temperature_celcius = temperature_celcius;
            self.pressure = pressure;
            self.composition = composition;
            self.which_pH_scale = pH_scale;
            self.which_ks = which_ks;
            self.co2_pressure_correction = co2_pressure_correction;

            self.pH_scale_conversion = [pHScale(pH_scale,composition,temperature_celcius,composition.salinity,which_ks,0.0),...
                                        pHScale(pH_scale,composition,temperature_celcius,composition.salinity,which_ks,pressure)];
        end
        function output = select(selected)
            output = KControls(self.temperature(selected),...
                               self.pressure(selected),...
                               self.composition.select(selected),...
                               self.which_pH_scale(selected),...
                               self.which_ks(selected),...
                               self.co2_pressure_correction(selected));
        end
    end
end