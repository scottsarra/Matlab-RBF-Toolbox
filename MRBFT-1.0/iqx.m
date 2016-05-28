%     Matlab Radial Basis Function Toolkit (MRBFT)
% 
%     Project homepage:    http://www.scottsarra.org/rbf/rbf.html
%     Contact e-mail:      sarra@marshall.edu
% 
%     Copyright (c) 2016 Scott A. Sarra
% 
%     Licensing: MRBFT is under the GNU General Public License ("GPL").
%     
%     GNU General Public License ("GPL") copyright permissions statement:
%     **************************************************************************
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


classdef iqx < rbfx
    methods
        function obj = iqx()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end
       
        function v = rbf(obj,r,s), v = 1./(1 + (s.*r).^2 ); end
       
        function d = D1(obj,r,s,x), d = -(2*x.*s.^2)./(1.0 + (s.*r).^2 ).^2; end
       
        function d = D2(obj, r, s, x)
            d = 2*s.^2.*(-r.^2.*s.^2 + 4*s.^2*x.^2 - 1)./(r.^2*s.^2 + 1).^3;
          %if nargin<5,  y=0;  end
	      %d = ( -2*s.^2 + s.^4.*(6*x.^2 - 2*y.^2)  )./(1.0 + (s.*r).^2 ).^3;
        end
        
        function d = D3(obj, r, s, x) 
            d = (-48*x.^3.*s.^6)./(r.^2.*s.^2 + 1).^4 + (24*x.*s.^4)./(r.^2.*s.^2 + 1).^3;
        end
        
        function d = D4(obj, r, s, x) 
            d = (384*x.^4.*s.^8)./(r.^2.*s.^2 + 1).^5 - (288*x.^2.*s.^6)./(r.^2.*s.^2 + 1).^4 + (24*s.^4)./(r.^2.*s.^2 + 1).^3; 
        end
        
        function d = G(obj, r, s, x, y)   % Gradient
           d = -2*s.^2.*(x + y)./(r.^2.*s.^2 + 1).^2;
        end
        
        function d = L(obj, r, s)         % Laplacian
           d = 4*s.^2.*(r.^2.*s.^2 - 1)./(1 + (s.*r).^2 ).^3;
        end
        
         % x and y not used but required by the abstract function definition in the superclass
        function d = B(obj, r, s, x, y)   % Biharmonic operator   
           d = 64*( s.^4 - 4*r.^2.*s.^6 + r.^4.*s.^8  )./(1 + (s.*r).^2 ).^5;
        end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y) 
           d = 8*x.*s.^4.*(1 - 5*y.^2.*s.^2 + x.^2.*s.^2  )./(1 + (s.*r).^2 ).^4;
       end
       
       function d = D22(obj, r, s, x, y) 
           d = -8*s.^4.*(-1 + 4*r.^2.*s.^2 + 5*(x.^4 + y.^4).*s.^4 - 38*x.^2.*y.^2.*s.^4  )./(1 + (s.*r).^2 ).^5;
       end
      
   end % methods
end  % class
