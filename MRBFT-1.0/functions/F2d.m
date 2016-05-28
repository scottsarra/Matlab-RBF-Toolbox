
% Franke function
% NOT FULLY IMPLEMENTED

classdef F2d < Function2d     

    methods(Static)
        function obj = F2d(), obj@Function2d(); end
        
        function f = F(x,y), f = 0.75.*exp(-0.25.*(9.*x-2).^2 - 0.25.*(9.*y-2).^2) + 0.75.*exp(-((9.*x+1).^2)./49 - ((9.*y+1).^2)./10) + ...
                    0.5.*exp(-0.25.*(9.*x-7).^2-0.25.*(9.*y-3).^2) -  0.2.*exp(-(9.*x-4).^2-(9.*y-7).^2);
        end
        
        function f = x1(x,y), f = 1;  end
        function f = x2(x,y), f = 1;  end
        function f = x3(x,y), f = 1;  end
        function f = x4(x,y), f = 1;  end
        
        function f = y1(x,y), f = 1;  end
        function f = y2(x,y), f = 1;  end
        function f = y3(x,y), f = 1;  end
        function f = y4(x,y), f = 1;  end
        
        function f = G(x,y), f = (9.0/1960)*(784*exp(-(4 - 9*x).^2 - (7-9*y).^2).*(-11 + 9*x + 9*y) - ...
                                490*exp(-(1.0/4)*(7 - 9*x).^2 - (9/4.0)*(1 - 3*y).^2).*(-10 + 9*x + 9*y) - ...
                                735*exp(-2 + 9*x - (81*x.^2)/4.0 + 9*y - (81*y.^2)/4.0).*(-4 + 9*x + 9*y) - ...
                                6*exp(-(1.0/49)*(1 + 9*x).^2 - (1.0/10)*(1 + 9*y).^2).*(59 + 90*x + 441*y));
        
        end
        function f = L(x,y), f = 1;  end
        
        function f = B(x,y), f = 1;  end
        
        function f = p12(x,y), f = 1;  end
        function f = p21(x,y), f = 1;  end
        function f = p22(x,y), f = 1;  end
        
    end % Static methods
    
end % class