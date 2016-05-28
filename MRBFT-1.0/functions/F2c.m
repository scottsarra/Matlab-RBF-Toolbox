classdef F2c < Function2d

    methods(Static)
        function obj = F2c(), obj@Function2d(); end
        
        function f = F(x,y), f = exp(x.*y); end
        
        function f = x1(x,y), f = exp(x.*y).*y;  end
        function f = x2(x,y), f = exp(x.*y).*y.^2;  end
        function f = x3(x,y), f = exp(x.*y).*y.^3;  end
        function f = x4(x,y), f = exp(x.*y).*y.^4;  end
        
        function f = y1(x,y), f = exp(x.*y).*x;  end
        function f = y2(x,y), f = exp(x.*y).*x.^2;  end
        function f = y3(x,y), f = exp(x.*y).*x.^3;  end
        function f = y4(x,y), f = exp(x.*y).*x.^4;  end
        
        function f = G(x,y), f = exp(x.*y).*(x + y);  end
        function f = L(x,y), f = exp(x.*y).*(x.^2 + y.^2);  end
        
        function f = B(x,y), f = exp(x.*y).*(4 + x.^4 + y.^4 + 8*x.*y + 2*x.^2.*y.^2);  end
        
        function f = p12(x,y), f = exp(x.*y).*x.*(2 + x.*y);  end
        function f = p21(x,y), f = exp(x.*y).*y.*(2 + x.*y);  end
        function f = p22(x,y), f = exp(x.*y).*(2 + 4*x.*y + x.^2.*y.^2);  end
        
    end % Static methods
    
end % class