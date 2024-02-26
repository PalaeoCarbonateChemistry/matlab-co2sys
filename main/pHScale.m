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
        function [pH_total,pH_seawater,pH_free,pH_NBS] = find_pH_on_all_scales(self,pH,Ks,selected,which_ks,salinity,temp_c)

            temp_k = temp_c+273.15;

            [K0,K1,K2,KW,KB,KF,KS,KP1,KP2,KP3,KSi,KNH4,KH2S] = Ks.unpack();
            fH = KFunctions.calculate_fH(which_ks,salinity,temp_k);
        
            fH_selected = fH(selected);
            FREEtoTOT = (1 + self.composition.sulphate(selected)./KS(selected)); % ' pH scale conversion factor
            SWStoTOT  = (1 + self.composition.sulphate(selected)./KS(selected))./(1 + self.composition.sulphate(selected)./KS(selected) + self.composition.fluorine(selected)./KF(selected));% ' pH scale conversion factor
            factor=nan(sum(selected),1);
            nF=self.which_pH_scale(selected)==1;  %'"pHtot"
            factor(nF) = 0;
            nF=self.which_pH_scale(selected)==2; % '"pHsws"
            factor(nF) = -log(SWStoTOT(nF))./log(0.1);
            nF=self.which_pH_scale(selected)==3; % '"pHfree"
            factor(nF) = -log(FREEtoTOT(nF))./log(0.1);
            nF=self.which_pH_scale(selected)==4;  %'"pHNBS"
            factor(nF) = -log(SWStoTOT(nF))./log(0.1) + log(fH_selected(nF))./log(0.1);
            pHtot  = pH    - factor;    % ' pH comes into this sub on the given scale
            pHNBS  = pHtot - log(SWStoTOT) ./log(0.1) + log(fH_selected)./log(0.1);
            pHfree = pHtot - log(FREEtoTOT)./log(0.1);
            pHsws  = pHtot - log(SWStoTOT) ./log(0.1);
            
            pH_total = pHtot;
            pH_seawater = pHsws;
            pH_free = pHfree;
            pH_NBS = pHNBS;
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