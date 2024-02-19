classdef Composition
    properties
        ammonia
        boron
        bromine
        calcium
        chloride
        fluorine
        magnesium
        potassium
        phosphate
        silicate
        sodium
        strontium
        sulphate
        sulphide

        salinity

        peng_correction
    end
    methods
        function self = Composition(salinity)
            self.salinity = salinity;

            self.ammonia = NaN(numel(self.salinity),1);
            self.boron = NaN(numel(self.salinity),1);
            self.bromine = NaN(numel(self.salinity),1);
            self.calcium = NaN(numel(self.salinity),1);
            self.chloride = NaN(numel(self.salinity),1);
            self.fluorine = NaN(numel(self.salinity),1);
            self.magnesium = NaN(numel(self.salinity),1);
            self.potassium = NaN(numel(self.salinity),1);
            self.phosphate = NaN(numel(self.salinity),1);
            self.silicate = NaN(numel(self.salinity),1);
            self.sodium = NaN(numel(self.salinity),1);
            self.strontium = NaN(numel(self.salinity),1);
            self.sulphate = NaN(numel(self.salinity),1);
            self.sulphide = NaN(numel(self.salinity),1);
        end

        function self = remove_freshwater_species(self,which_ks)
            selected = (which_ks.k1_k2==8);

            self.boron(selected) = 0;
            self.ammonia(selected) = 0;
            self.phosphate(selected) = 0;
            self.silicate(selected) = 0;
            self.sulphide(selected) = 0;
        end
        function self = adjust_geosecs_species(self,which_ks)
            selected = (which_ks.k1_k2==6);

            self.phosphate(selected) = 0;
            self.silicate(selected) = 0;
            self.sulphide(selected) = 0;
            self.ammonia(selected) = 0;

            selected = (which_ks.k1_k2==6 | which_ks.k1_k2==7);
            self.calcium(selected) = 0.01026.*self.salinity(selected)./35;
            self.boron(selected) = 0.0004106.*self.salinity(selected)./35;
        end
        function self = calculate_peng_correction(self,which_ks)
            self.peng_correction = zeros(numel(self.salinity),1); 

            selected = (which_ks.k1_k2==7); 
            self.peng_correction(selected) = self.phosphate(selected);
        end

        function self = set_ammonia_concentration(self,ammonia)
            self.ammonia = ammonia;
        end
        function self = set_boron_concentration(self,boron)
            self.boron = boron;
        end
        function self = set_bromine_concentration(self,bromine)
            self.bromine = bromine;
        end
        function self = set_calcium_concentration(self,calcium)
            self.calcium = calcium;
        end
        function self = set_chloride_concentration(self,chloride)
            self.chloride = chloride;
        end
        function self = set_fluorine_concentration(self,fluorine)
            self.fluorine = fluorine;
        end
        function self = set_magnesium_concentration(self,magnesium)
            self.magnesium = magnesium;
        end
        function self = set_potassium_concentration(self,potassium)
            self.potassium = potassium;
        end
        function self = set_phosphate_concentration(self,phosphate)
            self.phosphate = phosphate;
        end
        function self = set_silicate_concentration(self,silicate)
            self.silicate = silicate;
        end
        function self = set_sodium_concentration(self,sodium)
            self.sodium = sodium;
        end
        function self = set_strontium_concentration(self,strontium)
            self.strontium = strontium;
        end
        function self = set_sulphate_concentration(self,sulphate)
            self.sulphate = sulphate;
        end
        function self = set_sulphide_concentration(self,sulphide)
            self.sulphide = sulphide;
        end

        function boron = estimate_boron_from_salinity(self,which_boron)
            % selected=(which_ks.k1_k2~=6 & which_ks.k1_k2~=7 & which_ks.k1_k2~=8); % All other cases
            boron = self.boron;
            
            selected = isnan(self.boron) & which_boron==1;
            boron(selected) = 0.0004157.*self.salinity(selected)./35; % in mol/kg-SW
            
            selected = isnan(self.boron) & which_boron==2;
            boron(selected) = 0.0004326.*self.salinity(selected)./35; % in mol/kg-SW
        end
        function bromine = estimate_bromine_from_salinity(self)
            bromine = (0.0034730./79.904).*(self.salinity./1.80655);
        end
        function calcium = estimate_calcium_from_salinity(self)
            calcium = (0.02128./40.087).*(self.salinity./1.80655);
        end
        function chloride = estimate_chloride_from_salinity(self)
            chloride = (0.9989041./35.453).*(self.salinity./1.80655);
        end
        function fluorine = estimate_fluorine_from_salinity(self)
            fluorine = (0.000067./18.998).*(self.salinity./1.80655); % in mol/kg-SW
        end
        function magnesium = estimate_magnesium_from_salinity(self)
            magnesium =  (0.0662600./24.305).*(self.magnesium./1.80655);
        end
        function potassium = estimate_potassium_from_salinity(self)
            potassium = (0.0206000./39.0983).*(self.salinity./1.80655);
        end
        function sodium = estimate_sodium_from_salinity(self)
            sodium = (0.5564924./22.989769).*(self.salinity./1.80655);
        end
        function strontium = estimate_strontium_from_salinity(self)
            strontium = (0.0004100./87.62).*(self.salinity./1.80655);
        end
        function sulphate = estimate_sulphate_from_salinity(self)
            sulphate = (0.14./96.062).*(self.salinity./1.80655);
        end

        function self = estimate_all_from_salinity(self,which_boron)
            self.boron = self.estimate_boron_from_salinity(which_boron);
            self.calcium = self.estimate_calcium_from_salinity();
            self.chloride = self.estimate_chloride_from_salinity();
            self.fluorine = self.estimate_fluorine_from_salinity();
            self.magnesium = self.estimate_magnesium_from_salinity();
            self.potassium = self.estimate_potassium_from_salinity();
            self.sodium = self.estimate_sodium_from_salinity();
            self.strontium = self.estimate_strontium_from_salinity();
            self.sulphate = self.estimate_sulphate_from_salinity();
        end

        function [ammonia,boron,bromine,calcium,chloride,fluorine,magnesium,potassium,phosphate,silicate,sodium,strontium,sulphate,sulphide] = unpack(self)
            ammonia = self.ammonia;
            boron = self.boron;
            bromine = self.bromine;
            calcium = self.calcium;
            chloride = self.chloride;
            fluorine = self.fluorine;
            magnesium = self.magnesium;
            potassium = self.potassium;
            phosphate = self.phosphate;
            silicate = self.silicate;
            sodium = self.sodium;
            strontium = self.strontium;
            sulphate = self.sulphate;
            sulphide = self.sulphide;
        end
        function [ammonia,boron,fluorine,phosphate,silicate,sulphate,sulphide] = unpack_alkalinity(self)
            ammonia = self.ammonia;
            boron = self.boron;
            fluorine = self.fluorine;
            phosphate = self.phosphate;
            silicate = self.silicate;
            sulphate = self.sulphate;
            sulphide = self.sulphide;
        end
    end
end