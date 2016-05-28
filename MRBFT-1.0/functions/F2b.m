classdef F2b < Function2d

    methods(Static)
        function obj = F2b(), obj@Function2d(); end
        
        function f = F(x,y), f = exp(0.5*x + 0.2*y).*cos(x.*y); end
        
        function f = x1(x,y), f = 0.5*exp(0.5*x + 0.2*y).*( cos(x.*y) - 2*y.*sin(x.*y) );  end
        function f = x2(x,y), f = -0.25*exp(0.5*x + 0.2*y).*( (4*y.^2 - 1).*cos(x.*y) + 4*y.*sin(x.*y) );  end
        function f = x3(x,y), f = 0.125*exp(0.5*x + 0.2*y).*( (1 - 12*y.^2).*cos(x.*y) + 2*y.*(4*y.^2 - 3).*sin(x.*y) );  end
        function f = x4(x,y), f = 0.0625*exp(0.5*x + 0.2*y).*( (1 - 24*y.^2 + 16*y.^4).*cos(x.*y) + 8*y.*(-1 + 4*y.^2).*sin(x.*y) );  end
        
        function f = y1(x,y), f = 0.2*exp(0.5*x + 0.2*y).*( cos(x.*y) - 5*x.*sin(x.*y) );  end
        function f = y2(x,y), f = -0.04*exp(0.5*x + 0.2*y).*( (25*x.^2 - 1).*cos(x.*y) + 10*x.*sin(x.*y) );  end
        function f = y3(x,y), f = 0.008*exp(0.5*x + 0.2*y).*( (1 - 75*x.^2).*cos(x.*y) + 5*x.*(25*x.^2 - 3).*sin(x.*y) );  end
        function f = y4(x,y), f = 0.0016*exp(0.5*x + 0.2*y).*( (1 - 150*x.^2 + 625*y.^4).*cos(x.*y) + 20*x.*(-1 + 25*x.^2).*sin(x.*y) );  end
        
        function f = G(x,y), f = 0.1*exp(x/2 + y/5).*(7*cos(x.*y) - 10*(x + y).*sin(x.*y));  end
        function f = L(x,y), f = -0.01*exp(0.5*x + 0.2*y).*( (-29 + 100*x.^2 + 100*y.^2).*cos(x.*y) + 20*(2*x + 5*y).*sin(x.*y) );  end
        
        function f = B(x,y), f = 0.0001*exp(0.5*x + 0.2*y).*( (-39159 + 10000*x.^4 - 15800*y.^2 -16000*y + 10000*y.^4 - 8000*x.*(5+y) + 200*x.^2.*(100*y.^2 - 37)).*cos(x.*y) + ...
                                                             40*( 200*x.^3 + 500*x.^2.*y + x.*(-58 + 2000*y + 200*y.^2) + 5*(-40 - 29*y + 100*y.^3).*sin(x.*y) ) );  end
        
        function f = p12(x,y), f = -0.02*exp(0.5*x + 0.2*y).*( (-1 + 25*x.^2 + 20*x.*(5+y)).*cos(x.*y) + 2*(10 + 5*x + y - 25*x.^2.*y).*sin(x.*y) );  end
        function f = p21(x,y), f = -0.05*exp(0.5*x + 0.2*y).*( (-1 + 40*y + 20*x.*y + 4*y.^2).*cos(x.*y) + (20 + 5*x + 4*y - 20*x.*y.^2).*sin(x.*y) );  end
        function f = p22(x,y), f = 0.01*exp(0.5*x + 0.2*y).*( (-199 - 80*y - 4*y.^2 -200*x -40*x.*y - 25*x.^2 + 100*x.^2.*y.^2).*cos(x.*y) + 2*(50*x.^2.*y -2*(10+y) + 5*x.*(-1 + 40*y + 4*y.^2)).*sin(x.*y) );  end
        
    end % Static methods
    
end % class