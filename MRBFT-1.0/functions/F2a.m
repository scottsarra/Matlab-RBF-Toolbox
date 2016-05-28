classdef F2a < Function2d

    methods(Static)
        function obj = F2a(), obj@Function2d(); end
        
        function f = F(x,y), f = x.^3.*log(1+y) + y./(1+x); end
        
        function f = x1(x,y), f = 3*x.^2.*log(y + 1) - y./(x + 1).^2;  end
        function f = x2(x,y), f = 2*(3*x.*log(y + 1) + y./(x + 1).^3);  end
        function f = x3(x,y), f = 6*(-y./(x + 1).^4 + log(y + 1));  end
        function f = x4(x,y), f = 24*y./(x + 1).^5;  end
        
        function f = y1(x,y), f = (x.^3.*(x + 1) + y + 1)./((x + 1).*(y + 1));  end
        function f = y2(x,y), f = -x.^3./(y + 1).^2;  end
        function f = y3(x,y), f = 2*x.^3./(y + 1).^3;  end
        function f = y4(x,y), f = -6*x.^3./(y + 1).^4;  end
        
        function f = G(x,y), f = (x.^3.*(x + 1).^2 + 3*x.^2.*(x + 1).^2.*(y + 1).*log(y + 1) - y.*(y + 1) + (x + 1).*(y + 1))./((x + 1).^2.*(y + 1)) ;  end
        function f = L(x,y), f = -x.^3./(y + 1).^2 + 2*(3*x.*log(y + 1) + y./(x + 1).^3);  end
        function f = B(x,y), f = -6*x.^3./(y + 1).^4 - 12*x./(y + 1).^2 + 24*y./(x + 1).^5;  end
        function f = p12(x,y), f = -3*x.^2./(y + 1).^2;  end
        function f = p21(x,y), f = 2*(3*x./(y + 1) + (x + 1).^(-3));  end
        function f = p22(x,y), f = -6*x./(y + 1).^2;  end
        
    end % Static methods
    
end % class









