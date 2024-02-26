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
            
            fH = calculate_fH(which_ks,salinity,temp_k);

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