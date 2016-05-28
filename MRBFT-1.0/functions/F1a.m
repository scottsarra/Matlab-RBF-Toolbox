classdef F1a < Function1d

    methods(Static)
        function obj = F1a(), obj@Function1d(); end
        
        function f = F(x), f = exp(sin(pi*x)); end    % [-1,1]
        
        function f = x1(x)
            f = pi*cos(pi*x).*exp(sin(pi*x));
        end
        
        function f = x2(x)
            f =pi*pi*( -sin(pi*x) + cos(pi*x).^2 ).*exp(sin(pi*x));
        end
        
        function f = x3(x)
            f = -pi^3*(sin(pi*x) + 3).*exp(sin(pi*x)).*sin(pi*x).*cos(pi*x); 
        end
        
        function f = x4(x)
            f = pi^4*(sin(pi*x).^4 + 6*sin(pi*x).^3 + 5*sin(pi*x).^2 - 5*sin(pi*x) - 3).*exp(sin(pi*x)); 
        end
        
    end % Static methods
    
end % class
   