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


classdef rbfx

% ----------------------------------------------------------------------------
% ----------------------- Abstract methods------------------------------------
%             used to define a common interface for all subclasses
%                    must be implemented by all subclasses.
% ----------------------------------------------------------------------------

   methods(Abstract = true)    
       v = rbf(obj,r,s);         % RBF definition
       d = D1(obj,r,s,x);        % first derivative wrt x
       d = D2(obj, r, s, x);     % second derivative wrt x
       d = D3(obj, r, s, x);     % third derivative wrt x
       d = D4(obj, r, s, x);     % fourth derivative wrt x
       d = G(obj, r, s, x, y);   % Gradient
       d = L(obj, r, s);         % Laplacian
       d = B(obj, r, s, x, y);   % Biharmonic operator
       d = D12(obj, r, s, x, y); % mixed partial derivative 
       d = D22(obj, r, s, x, y); % mixed partial derivative        
   end
   
% ----------------------------------------------------------------------------
% ---------------- static methods --------------------------------------------
% ----------------------------------------------------------------------------
    
    methods(Static)
         
% ----------------------------------------------------------------------------
% ------------------ distance matrices ---------------------------------------
% ----------------------------------------------------------------------------

% distanceMatrix1d
% inputs
%   xc   N x 1 vector of centers
%   x    M x 1 vector of evaluation points (optional)
%
% outputs
%   r    signed distance matrix
%        N x N,   r_{ij} = dist between center i and j if called as distanceMatrix1d(xc)
%        M x N,   r_{ij} = dist between evaluation point i and center j if called as distanceMatrix1d(xc,x)
%        N x N/2  if called as r = phi.distanceMatrix1d(xc(1:N/2),xc)
%                 returns the left half of the distance matrix needed
%                 for a centro center distribution
%
% example usage:
%       1) condVaccuracy.m, 2) rbfInterpConvergenceB.m
         
    function r = distanceMatrix1d(xc,x)
        xc = xc(:);                   %  make sure xc is a column vector
        o = ones(1,length(xc));
        if nargin==1
            r = xc*o;
            r = r - r';  
        else
            x = x(:);
            r = x*o - ones(length(x),1)*xc';
        end
    end
     

% distanceMatrix2d
%
% inputs
%   xc   N x 1 vectors of centers  XC = (xc,yc)
%   yc
%   x    M x 1 vectors of evaluation points X = (x,y) (optional)
%   y
%
% outputs
%   r    signed distance matrix
%        N x N,   r_{ij} = dist between center i and j if called as distanceMatrix1d(xc,yc)
%        M x N,   r_{ij} = dist between evaluation point i and center j if called as distanceMatrix1d(xc,yc,x,y)
%        N x N/2  if called as r = phi.distanceMatrix1d(xc(1:N/2),yc(1:N/2),xc,yc)
%                 returns the left half of the distance matrix needed
%                 for a centro center distribution
%
% example usage:
%     1) diffusionReactionCentro.m, 2) mdiRegularization.m, 3) poissonCentro.m
         
       
    function [r, rx, ry] = distanceMatrix2d(xc,yc,x,y)
        xc = xc(:);   yc = yc(:);
        o = ones(1,length(xc));
        if nargin==2
            rx = (xc*o - (xc*o)');
            ry = (yc*o - (yc*o)');
            r = sqrt( rx.^2  + ry.^2 ); 
        else
            om = ones(length(x),1);
            x = x(:);   y = y(:);
            rx = (x*o - om*xc');
            ry = (y*o - om*yc');
            r = sqrt( rx.^2 + ry.^2 ); 
         end
    end
       
 
% distanceMatrix3d
%
% inputs
%   xc   N x 1 vectors of centers  XC = (xc,yc,zc)
%   yc
%   zc
%   x    M x 1 vectors of evaluation points X = (x,y,z) (optional)
%   y
%   z
%
% outputs
%   r    signed distance matrix
%        N x N,   r_{ij} = dist between center i and j if called as distanceMatrix1d(xc,yc,zc)
%        M x N,   r_{ij} = dist between evaluation point i and center j if called as distanceMatrix1d(xc,yc,zc,x,y,z)
%        N x N/2  if called as r = phi.distanceMatrix1d(xc(1:N/2),yc(1:N/2),zc(1:N/2),xc,yc,zc)
%                 returns the left half of the distance matrix needed
%                 for a centro center distribution
%
% example usage: 1) interp3d.m, 2) interp3dCentro.m
       
    function [r, rx, ry, rz] = distanceMatrix3d(xc,yc,zc,x,y,z)
        xc = xc(:);   yc = yc(:);  zc = zc(:);
        o = ones(1,length(xc));
        if nargin==3
            rx = (xc*o - (xc*o)');
            ry = (yc*o - (yc*o)');
            rz = (zc*o - (zc*o)');
            r = sqrt( rx.^2  + ry.^2 + rz.^2 ); 
        else
            om = ones(length(x),1);
            x = x(:);   y = y(:);  z = z(:);
            rx = (x*o - om*xc');
            ry = (y*o - om*yc');
            rz = (z*o - om*zc');
            r = sqrt( rx.^2 + ry.^2 + rz.^2 ); 
         end
    end
  
% ----------------------------------------------------------------------------
% ---------- regularized SPD linear system solvers ---------------------------
% ----------------------------------------------------------------------------
       
% solve - solves the SPD linear system B a = f for a with the option to regularize
%         by the method of diagonal increments (MDI)
% inputs 
%    B    N x N symmetric positive definite (SPD) matrix
%    f    N x 1 vector
%   mu    (optional) MDI regularization parameter.  Use mu = 0 for no regularization
%  safe   true - uses backslash with error checking etc.
%         false - uses a Cholesky factorization.  Faster, but it the matrix is 
%                 severely ill-conditioned and/or the regularization parameter is too
%                 small the matrix may fail to be numerically SPD and the Cholesky
%                 factorization will fail
%
% outputs
%    a    N x 1 solution vector
%
% example usage:
%   1) mdiExample.m, 2) mdiRegularization.m, 3) rbfInterpConvergence.m

    function a = solve(B,f,mu,safe)  
        if ~exist('mu','var'), mu = 5e-15;  end
        if ~exist('safe','var'), safe = true;  end
  
        if mu>0
            N = length(f);
            B(1:N+1:end) = B(1:N+1:end) + mu;  % C = B + mu*eye(N);
        end
          
        if safe
            a = B\f;
        else
            L = chol(B,'lower');
            a = L'\( L\f );      
       end
    end
            
% -------------------------------------------------------------------------       
 
% dm - forms the deriviative matrix D = H*inv(B) by solving the system D B = H for D
% inputs
%    B   N x N SPD system matrix
%    H   N x N derivative evaluation matrix
%   mu    (optional) MDI regularization parameter.  Use mu = 0 for no regularization
%  safe   true - uses backslash with error checking etc.
%         false - uses a Cholesky factorization.  Faster, but if the matrix is 
%                 severely ill-conditioned and/or the regularization parameter is too
%                 small the matrix may fail to be numerically SPD and the Cholesky
%                 factorization will fail
%
% outputs
%    D   N x N differentiation matrix
%
% example usage:
%  1) diffusionReactionCentro.m

    function D = dm(B,H,mu,safe)      
        if ~exist('mu','var'), mu = 5e-15;  end
        if ~exist('safe','var'), safe = true;  end
        
        if mu>0
            s = size(B); N = s(1);
            B(1:N+1:end) = B(1:N+1:end) + mu;  % B = B + mu*eye(N);
        end
          
        if safe
            D = H/B; 
        else 
            L = chol(B,'lower');
            D = (L'\(L\H'))';   
       end
    end
       

% ----------------------------------------------------------------------------   
% -------------- variable shape parameters -----------------------------------
% ----------------------------------------------------------------------------

% variableShape 
%
% example usage: variableShapeInterp1d.m
%
%  inputs
%   sMin    minimum value of the shape parameter
%   sMax    maximum value of the shape parameter
%    N      number of columns
%    M      number of rows
%  opt  1, exponentially varying (Kansa)
%          Computers and Mathematics with Applications v. 24, no. 12, 1992. 
%       2, linearly varying
%       3, randonly varying (Sarra and Sturgil)
%          Engineering Analysis with Boundary Elements, v. 33, p. 1239-1245, 2009.
%
%  outputs
%    s1   N x N matrix with constant shapes in each column
%         call as, s1 = rbfx.variableShape(sMin,sMax,N)
%    s2   M x N matrix with constant shapes in each column 
%         (optional, for interpolation evaluation matrix) 
%          call as, [s1, s2] = rbfx.variableShape(sMin,sMax,N,M)
%
% example usage:
%   1) variableShapeInterp1d.m

    function varargout = variableShape(sMin,sMax,opt,N,M)
        if nargin<5, M = []; end
        nOutputs = nargout;
        varargout = cell(1,nOutputs);

        if opt==1          
            sMin = sMin^2;      sMax = sMax^2;
            s = sqrt( sMin*(sMax/sMin).^((0:N-1)./(N-1)) );      
        elseif opt==2   
            s = sMin + ((sMax - sMin)/(N-1)).*(0:N-1);
        else     
            s = rand(1,N);
            s = sMin + (sMax - sMin)*s; 
        end
        
        if nOutputs==2
            varargout{1} = repmat(s,N,1);  
            varargout{2} = repmat(s,M,1);
        else
           varargout{1} = repmat(s,N,1);  
        end
    end 



% -------------------------------------------------------------------------


end % methods


    
% --------------------------------------------------------------------------- 
    
   
end % classdef
