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

classdef rbfCentro


   
% ----------------------------------------------------------------------------
% ---------------- static methods --------------------------------------------
% ----------------------------------------------------------------------------
    
    methods(Static)
         
% ----------------------------------------------------------------------------
% ------------ centro symmetric functions ------------------------------------
% ----------------------------------------------------------------------------


% Note: In testing a matrix for centro symmetry, the full N x N matrix constructed
%       by non-centrosymmetric algorithms must be used as extending an arbitrary 
%       N x (N/2) matrix to an N x N matrix with rbfCentro.fullCentroMatrix will make
%       the full matrix centrosymmetric


% isCentro and isSkewCentro - test a matrix for  (skew)centrosymmetry
%
% input
%    B    a N x N matrix
%
% output
%    s    returns possibly zero or more probably a very small number, e.g. 1e-15, 
%         for a N x N centrosymmetric input. otherwise, the return of a larger
%         number indicates the input matrix is
%         not centrosymmetric


    function c = isCentro(B)   
        c = max( max( abs( B  - fliplr(flipud(B))  ) ));
    end


   function c = isSkewCentro(B)  
       c = max( max( abs( B  + fliplr(flipud(B))  ) ));
     %   c = max( max( abs(B  + rot90(B,2) ) ) );    % works also
   end

% hasSymmetry - tests a N x N matrix for both centrosymmetry and skew-centrosymmetry
%
% input
%    B    a N x N matrix
%
% output
%    text message written to the Matlab command window
%
% example usage:
%   1) isCentroTest.m

   function hasSymmetry(B)
      c1 = rbfCentro.isCentro(B);
      c2 = rbfCentro.isSkewCentro(B);
      if c1 < 100*eps
          disp('centrosymmeric')
      end
      if c2 < 100*eps
          disp('skew centrosymmeric')
      end
      
      if (c1 > 100*eps) && (c2 > 100*eps)
       fprintf('no centro symmetry (c = %4.1e, s = %4.1e)\n',c1,c2);
      end
      
   end
   
   
   
   
 
    
%  fullCentroMatrix
%
%  Corrects a full centro matrix constructed with a standard algorithm so that is has
%  the correct symmetry or takes the left half of a matrix and expands it to have symmetry
%
%  inputs
%   Dh      either a N x N  centro or skew-centro symmetric matrix
%             or the or the (N) x (N/2) left half of a such a matrix
%  skew    1 (or any positve odd integer) if skew-centro
%          2 (or any positive even integer) if centro
%    N     the number of rows of Dh (must be passed as a mp object for extended precision to work)
%          NOTE: N must be even
%
% outputs
%    D     the full N x N matrix
%  
    
    function D = fullCentroMatrix(Dh,N,skew)     
        n = N/2;     
        D = zeros(N,N);                                         
        D(:,1:n) = Dh(:,1:n);                                           % left half
        D(1:n,n+1:N) = (-1)^(skew)*fliplr( flipud(  D( n+1:N, 1:n) ) );   %   D12
        D(n+1:N, n+1:N) = (-1)^(skew)*fliplr( flipud(  D(1:n, 1:n) ) );   %   D22
    end
    
% -------------------------------------------------------------------------

% solveCentro.  Solves a centrosymmetric linear system Ba = f.
%
%  inputs
%  B    either a N x N  centrosymmetric matrix (N must be even)
%       or the or the (N) x (N/2) left half of a centrosymmetric matrix
%  f    N x 1 vector
%  mu   MDI regularization parameter (mu=0 for no regularization)
% safe  true -> use backslash operator
%       false -> directly use Cholesky factorization
%
% outputs
%  a    N x 1 solution of B a = f
%
% example usage: 1) poissonCentro.m
              
    function a = solveCentro(B,f,mu,safe)   
        if ~exist('mu','var'), mu = 5e-15;  end
        if ~exist('safe','var'), safe = true;  end

        N = length(f);
        n = int64( N/2 );     
          
        A = B(1:n,1:n);                    % B11
        t1 = flipud( B(n+1:N,1:n) );       % J*B21
        L = A - t1;                        % B11 - J*B21
        M = A + t1;                        % B11 + J*B21
          
        if mu>0                                                          
            L(1:n+1:end) = L(1:n+1:end) + mu;  % L = L + mu*eye(n);
            M(1:n+1:end) = M(1:n+1:end) + mu;
        end
        
        b1 = f(1:n);
        t2 = flipud( f(n+1:N) );        %  t2 = J*b2;
        b1p = b1 - t2;
        b2p = b1 + t2;
          
        if safe
           x1h = L\b1p;
           x2h = M\b2p;
        else
            L1 = chol(L,'lower');
            x1h = L1'\(L1\b1p);
   
            L2 = chol(M,'lower');
            x2h = L2'\(L2\b2p);
        end

        a = vertcat( 0.5*( x1h + x2h ), 0.5*flipud(x2h - x1h));         
    end
       

% -------------------------------------------------------------------------

% centroDM.  Constructs a RBF differentiation matrix.
%
% inputs
%    B    N x (N/2) left half (or full N x N) centro system matrix.
%         Only the left half is used.
%    F    N x (N/2) left half (or full N x N) derivative evaluation matrix.
%         Only the left half is used.
%    N    Number of columns in B and F. Must be passed as mp('N') for 
%         extended precision calculations.  N must be even.
%  rho    derivative order. 
%   mu    MDI regularization parameter (optional).  Should be passed for 
%         extended precisions as the default is for double precision.
%  safe   true (backslash), false (Cholesky)
%
% outputs
%    D    N x (N/2) left half of the DM. If the full DM is needed it can be
%         constructed with rbfCentro.fullCentroMatrix.
%
% example usage:
%    1) diffusionReactionCentro.m

    function D = centroDM(B,F,N,rho,mu,safe)   
        if ~exist('mu','var'), mu = 5e-15;  end
        if ~exist('safe','var'), safe = true;  end

        F = F*(-1)^(rho);  
        
        n = N/2;
        D = zeros(N,n);

        A = B(1:n,1:n);               % B11
        t1 = flipud( B(n+1:N,1:n) );  % J*B21
          
        L = A - t1;
        M = A + t1;
          
        if mu>0
            L(1:n+1:end) = L(1:n+1:end) + mu;  % L = L + mu*eye(n);
            M(1:n+1:end) = M(1:n+1:end) + mu;
        end
          
        b1 = F(1:n,1:n);
        t2 = flipud( F(n+1:N,1:n) );
        b1p = b1 - t2;
        b2p = b1 + t2;
          
        if safe
            x1h = L\b1p;
            x2h = M\b2p;
        else
            L1 = chol(L,'lower');
            L2 = chol(M,'lower');
          
            x1h = L1'\( L1\b1p );
            x2h = L2'\( L2\b2p );
        end
        
        D(1:n,1:n) = 0.5*( x1h + x2h )';                                        % D11
        D(n+1:N,1:n) = (-1)^rho*fliplr( flipud(  0.5*flipud(x2h - x1h)' ) );    % D21
    end
    
% ------ Condition number of a centrosymmetric system matrix ---------

% inputs
%    B    N x (N/2) left half (or full N x N) centro system matrix.
%         Only the left half is used.
%   mu    MDI regularization parameter
%
% outputs:
%    kappaL
%    kappaM
%    kappaB   2-norm condition number of the matrix B
%
% example usage:
%     1) diffusionReactionCentro.m, 2) poissonCentro.m

    function [kappaB, kappaL, kappaM] = centroConditionNumber(B,mu)
       [L,M] = rbfCentro.centroDecomposeMatrix(B,0);
       N = size(L);  
       I = eye(N(2));
       sL = svd(2*(L + mu*I));  sM = svd(2*(M + mu*I));
       kappaL = max(sL)/min(sL);     kappaM = max(sM)/min(sM);
       s = [sL; sM];
       kappaB = max(s)/min(s);
    end
    
% NOTE: Mathematically, the following funtion that uses eigenvalues rather
%       than singular values is equivalent.  However, if B is very ill-conitioned
%       for example cond(B)>10e15, the function using eig may return a complex
%       number as a condition number whereas the one using the SVD will always
%       return a real number.


    function kappa = centroConditionNumberEig(B,mu)
       [L,M] = rbfCentro.centroDecomposeMatrix(B,0);
       N = size(L);  I = eye(N(2));
       ew = [eig(2*(L + mu*I)); eig(2*(M + mu*I))];
       kappa = max(ew)/min(ew);
    end


% ---------- Parity Matrix Multiplication (N even) ------------------------

% centroDecomposeMatrix.  A (skew) centrosymmetric matrix is similar to a 
%                         block diagonal matric
%                          [ L   O ]
%                          [ 0   M ]
%  Given D, this functions computes L and M which are used in condition
%  number, eigenvalue, and multiplication algorithsm.
%
%  inputs
%  D    either a N x N  centrosymmetric matrix DM
%       or the or the (N) x (N/2) left half of a centrosymmetric matrix DM
%
%  outputs
%    L/2    the even (N/2) x (N/2) DM
%    M/2    the odd  (N/2) x (N/2) DM
%
%  example usage:
%     1) diffusionReactionCentro.m, 2) poissonCentro.m

    function [L,M] = centroDecomposeMatrix(D,rho)      
        s = size(D);  
        N = s(1);     % number of rows
        N2 = N/2;     % number of cols
        
        a = D(1:N2,1:N2);                               % B11 
        b = (-1)^(rho+1)*flipud( D(N2+1:N,1:N2) );      % J*B21
        
        L = 0.5*( a - b);  % L = 0.5*(B11 - J*B21) or Lh = 0.5*( B11 + J*B21)
        M = 0.5*( a + b);  % M = 0.5*(B11 + J*B21) or Mh = 0.5*( B11 - J*B21 )
     
    end

% centroMult. Matrix-vector multiplication with a (skew)centrosymmetric matrix
%  
% inputs
%   u    N x 1 vector
%   L    the even (N/2) x (N/2) DM  (L and M from rbfCentro.centroDecomposeMatrix)
%   M    the odd  (N/2) x (N/2) DM
%
% outputs
%  ua    N x 1 vector that is the results of D*f
%
% example usage:
%    1) diffusionReactionCentro.m, 2) poissonCentro.m
    
    function ua = centroMult(u,L,M,rho)   % rho even, centro; rho odd, skew-centro
        u = u(:);  N = length(u);
        N2 = N/2;  k = 1:N2;
        
        t = flipud( u(N2+1:end) );     % decompose
        e = u(1:N2) + t;                  % xe_1           
        o = u(1:N2) - t;                  % xo_1
        
        uae = L*e;                       % fe_1
        uao = M*o;                       % fo_1
        
        ua = zeros(N,1);              % reconstruct
        ua(1:N2) = uae + uao;                      % f_1 = fe_1 + fo_1
        s1 = (-1)^(rho);  s2 = (-1)^(rho+1);
        ua(N2+1:end) = flipud( s1*uae + s2*uao );  % f_2 = s1*J*fe_1 + s2*J*fo_1

         
    end


    
% -------------------------------------------------------------------------
% ------- centrosymmetric center distributions in 2d domains --------------
% -------------------------------------------------------------------------

% inputs
%    x, y      centers covering an entire domain
%  symType     0 - y axis
%              1 - x axis
%              2 - origin
%    plt       logical variable to plot centers
%
% outputs      
%   xc, yc      centrosymmetric center distribution
%
% example usage:
%   1) diffusionReactionCentroDriver.m

    function [xc,yc] = centroCenters(x,y,symType,plt)
         x = x(:);  y = y(:);  % ensure column vectors
        
        if symType == 0        
            I = find(x>0);
            xc = [x(I); flipud(-x(I))];       
            yc = [y(I); flipud(y(I))];
        elseif symType == 1    
            I = find(y>0);
            xc = [x(I); flipud(x(I))];
            yc = [y(I); flipud(-y(I))];
        elseif symType == 2       
            I = find(y>x);
            xc = [x(I); flipud(-x(I))];
            yc = [y(I); flipud(-y(I))];
        end
        
        if plt, scatter( xc, yc,'b.'); end
       
    end


   
% -------------- symmetric center distributions on a circle ---------------

% centroCircle - centroCenters for the specific cace of a circle. Quasirandom 
%                centers on a circle of radius R with centro placement. The centers
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
%    xc, yc     center coordinates
%
%  example usage:
%   1) isCentroTest.m

    function [xc,yc] = centroCircle(N,cluster,ch,R,plt)
        if ~exist('plt','var'), plt = false;  end
        if ~exist('R','var'), R = 1;  end
        [x, y] = rbfCenters.circleCenters(N,cluster,ch,R,plt);
        x = x(:);  y = y(:);  
        
        I = find(y>x);        % extend about the origin
        xc = [x(I); flipud(-x(I))];
        yc = [y(I); flipud(-y(I))];
        
        if plt, scatter( xc, yc,'b.'); end
    end


end % methods


    
% --------------------------------------------------------------------------- 
    
   
end % classdef
