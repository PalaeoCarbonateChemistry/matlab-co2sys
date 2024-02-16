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
    end
end