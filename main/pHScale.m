classdef pHScale
    properties
        which_pH_scale
        composition

        ks
        kf

        seawater_to_total
        free_to_total
    end
    methods
        function self = pHScale(which_pH_scale,composition,temp_c,salinity,which_ks,pressure_bar)
            self.which_pH_scale = which_pH_scale;
            self.composition = composition;
            
            self.ks = KFunctions.calculate_ks(temp_c,pressure_bar,salinity,which_ks);
            self.kf = KFunctions.calculate_kf(temp_c,pressure_bar,salinity,which_ks);

            self.seawater_to_total = pHScale.calculate_seawater_to_total(composition,self.ks,self.kf);
            self.free_to_total = pHScale.calculate_free_to_total(composition,self.ks);
        end
        function pH_factor = calculate_pH_factor(self,temp_c,salinity,which_ks)
            temp_k = temp_c+273.15;
            
            fH = KFunctions.calculate_fH(which_ks,salinity,temp_k);

            % FindpHScaleConversionFactor:
            % this is the scale they will be put on
            pH_factor = NaN(numel(temp_c),1);

            selected = (self.which_pH_scale==1); % Total
            pH_factor(selected) = self.seawater_to_total(selected);

            selected = (self.which_pH_scale==2); % SWS, they are all on this now
            pH_factor(selected) = 1;

            selected = (self.which_pH_scale==3); % pH free
            pH_factor(selected) = self.seawater_to_total(selected)./self.free_to_total(selected);

            selected = (self.which_pH_scale==4); % pH NBS
            pH_factor(selected) = fH(selected);
        end

        function [pH_total,pH_seawater,pH_free,pH_NBS] = find_pH_on_all_scales(self,pH,K_controls)

            temp_k = K_controls.temperature_celcius + 273.15;

            fH = KFunctions.calculate_fH(K_controls.which_ks,K_controls.composition.salinity,temp_k);
        
            fH_selected = fH;

            pH_factor = NaN(numel(pH),1);

            current_pH_scale = self.which_pH_scale==1;  %'"pHtot"
            pH_factor(current_pH_scale) = 0;

            current_pH_scale=self.which_pH_scale==2; % '"pHsws"
            pH_factor(current_pH_scale) = -log(self.seawater_to_total(current_pH_scale))./log(0.1);

            current_pH_scale=self.which_pH_scale==3; % '"pHfree"
            pH_factor(current_pH_scale) = -log(self.free_to_total(current_pH_scale))./log(0.1);

            current_pH_scale=self.which_pH_scale==4;  %'"pHNBS"
            pH_factor(current_pH_scale) = -log(self.seawater_to_total(current_pH_scale))./log(0.1) + log(fH_selected(current_pH_scale))./log(0.1);
            
            pH_total  = pH    - pH_factor;    % ' pH comes into this sub on the given scale
            pH_NBS  = pH_total - log(self.seawater_to_total) ./log(0.1) + log(fH_selected)./log(0.1);
            pH_free = pH_total - log(self.free_to_total)./log(0.1);
            pH_seawater  = pH_total - log(self.seawater_to_total) ./log(0.1);            
        end
    end
    methods (Static=true)
        function seawater_to_total = calculate_seawater_to_total(composition,ks,kf)
            seawater_to_total = (1 + composition.sulphate./ks)./(1 + composition.sulphate./ks + composition.fluorine./kf);
        end
        function free_to_total = calculate_free_to_total(composition,ks)
            free_to_total =  1 + composition.sulphate./ks;
        end
    end
end