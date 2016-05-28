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

classdef rbfCenters

    
% ----------------------------------------------------------------------------
% ---------------- static methods --------------------------------------------
% ----------------------------------------------------------------------------
    
    methods(Static)
         
% Hammersley2d. Quasi-random Hammersley points on the unit square.
%
% inputs
%   N      number of centers
%  plt     logical variable: true -> plot centers, false -> no plot
%
% output
%   (x,y)  center locations
%
% example usage:
%   1) centroCenters.m

    function [x, y] = Hammersley2d(N,plt)
        if ~exist('plt','var'), plt = false;  end
        x = zeros(N,1);  y = zeros(N,1);
        
        for k = 0:N-1
            u = 0;
            p = 0.5;
            kk = k;
            while kk>0
                if bitand(kk,1)
                    u = u + p;
                end
                p = 0.5*p;
                kk = bitshift(kk,-1);
            end
            v = (k + 0.5)/N;
            x(k+1) = u;
            y(k+1) = v;
        end
        
        if plt, scatter(x,y,'b.'); end
        
    end


% Halton2d. Quasi-random Halton points on the unit square.
%
% inputs
%   N      number of centers
%  plt     logical variable: true -> plot centers, false -> no plot
%
% output
%   (x,y)  center locations

    function [x, y] = Halton2d(N,plt)
        if ~exist('plt','var'), plt = false;  end
        x = zeros(N,1);y = zeros(N,1);
        
        k = 0;
        while (k+1) ~= N
            u = 0;
            p = 0.5;
            kk = k;
            
            while kk>0
                if bitand(kk,1)
                    u = u + p;
                end
                p = 0.5*p;
                kk = bitshift(kk,-1);
            end
            
            v = 0;
            p2 = 3.0;        % prime2 which is taken to be 3
            ip = 1.0/p2;           
            p = ip;
            kk = k;
            
            while kk>0
                a = rem(kk,p2);
                if a~=0, v = v + a*p; end
                p = p*ip;
                kk = floor(kk/p2);
            end
            
            x(k+1) = u;
            y(k+1) = v;
            k = k + 1;
        end
    
        if plt, scatter(x,y,'b.'); end
    end

% squareCenters
%
% Quasirandom centers on a square [a,b]^2. The centers
% are either based on a Hammersley or Halton sequence.
%
%   inputs
%      N      Number of centers in the covering square.  The number of centers
%             returned is less than N.
%    cluster  Logical variable for clustering option
%     ch          1        Halton
%             otherswise   Hammersley
%     a, b    square [a,b] x [a,b]
%     plt     Logical variable for plotting option
%
%  outputs
%    x, y     center coordinates

function [x, y] = squareCenters(N,a,b,cluster,ch,plt)  
    if ~exist('plt','var'), plt = false;  end
    
    if ch == 1
        [x, y] = rbfCenters.Halton2d(N,false);
    else
        [x, y] = rbfCenters.Hammersley2d(N,false);
    end
    
    x = 2*x - 1;             % [0,1]^2 --> [-1,1]^2
    y = 2*y - 1;
    
    if cluster
        x = sin(0.5*pi*x);      % cluster
        y = sin(0.5*pi*y);
    end
    
    x = 0.5*(b-a)*x + 0.5*(b + a);        % [-1,1]^2 --> [a,b]^2
    y = 0.5*(b-a)*y + 0.5*(b + a);
    
    if plt, scatter(x,y,'b.'); end
    
end
    
    

    
   
% circleCenters - Quasirandom centers on a circle of radius R. The centers
%                 are either based on a Hammersley or Halton sequence.
%   inputs
%      N      Number of centers in the covering square.  The number of centers
%             returned is less than N.
%    cluster  Logical variable for clustering option
%     ch          1        Halton
%             otherswise   Hammersley
%      R      Radius of the circle
%     plt     Logical variable for plotting option
%
%  outputs
%    x, y     center coordinates
%
%  example usage: 1) interp2d_d.m
    
    
    function [x, y] = circleCenters(N,cluster,ch,R,plt)
        if ~exist('plt','var'), plt = false;  end
        if ~exist('R','var'), R = 1;  end
        if ch == 1
            [x, y] = rbfCenters.Halton2d(N);
        else
            [x, y] = rbfCenters.Hammersley2d(N);
        end
            
        x = 2*x - 1;                      % [0,1]^2 --> [-1,1]^2
        y = 2*y - 1;                                 
        I = find( x.^2 + y.^2 <= 1 );     % restrict from square to circle
        x = x(I);   y = y(I);
        
        if cluster
            [t,r] = cart2pol(x,y);
            r = sin(0.5*pi*r);
            [x,y] = pol2cart(t,r);
        end
        
        x = R*x;  y = R*y;                 % adjust to have radius R
        if plt, scatter(x, y,'b.'); end
    end
    
    
    
    
    
% circleUniformCenters - Uniform centers on a circle of radius R.
%
%    inputs 
%      N      The number of centers returned
%      R      Radius of the circle
%     plt     Logical variable for plotting option
%
%  outputs
%    x, y     center coordinates
%      Nb     The number of center located on the boundary which are in
%             the last Nb locations of the returned vector.  Useful for
%             enforcing PDE boundary conditions. To plot boundary centers: 
%                   hold on
%                   scatter( x(end-Nb+1:end), y(end-Nb+1:end),'ro')
%
% example usage:
%     1) poissonCentro.m

    function [x,y, Nb] = circleUniformCenters(N,R,plt)
        if ~exist('plt','var'), plt = false;  end
        if ~exist('R','var'), R = 1;  end
        x(1) = 0; y(1) = 0;
        Ns = round( (sqrt(pi+4*(N-1)) - sqrt(pi)) /(2*sqrt(pi)) );       % number of circles
        K = pi*(Ns+1)/(N-1);            % constant used in each loop

        for i = 1:Ns-1
            ri = i/Ns;
            ni = round(2*pi*ri/K);      % number of points on circle i
            t = linspace(0, 2*pi, ni+1)';   t = t(1:ni);

            if mod(i,2)==0,   % stagger the start of every other circle
                dt = t(2) - t(1);
                t = t + 0.5*dt;
            end

            x = [x; ri*cos(t)];   y = [y; ri*sin(t)];  
        end

        Nb = N - length(x);   % remaining points to be placed on the outter circle
        t = linspace(0, 2*pi, Nb + 1)';  t = t(1:Nb);
        x = [x; cos(t)];  y = [y; sin(t)];
        x = R*x;  y = R*y;                 % adjust to have radius R

        if plt, scatter( x, y,'b.'); end
    end
    

% -------------------------------------------------------------------------


end % methods

    
% --------------------------------------------------------------------------- 
    
   
end % classdef
