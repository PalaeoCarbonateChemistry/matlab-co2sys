classdef WhichKs
    properties
        k1_k2
        kso4
        kf
    end
    methods
        function self = WhichKs(k1_k2,kso4,kf)
            self.k1_k2 = k1_k2;
            self.kso4 = kso4;
            self.kf = kf;
        end
        function output = select(self,selected)
            if any(selected)
                output = WhichKs(self.k1_k2(selected),...
                                 self.kso4(selected),...
                                 self.kf(selected));
            else
                output = NaN;
            end
        end
    end
end