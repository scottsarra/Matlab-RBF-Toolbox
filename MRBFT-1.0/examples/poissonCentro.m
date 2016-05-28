% poissonCentro.m
%
% Solves a 2d steady PDE problem, a Poisson equation, on a circular domain 
% with Dirhclet boundary conditions, u_xx + u_yy = f(x,y).
%
% The problem is solved two ways: 1) standard algorithms, 2) centrosymmetric 
% algorithms.  With N = 5000 the accuracy of the 2 approaches is the same
% but the centrosymmetric approach is approximately 5 times faster and 
% requires half the storage.


clear, home, close all

CENTRO = false
mu = 0;
phi = iqx();

%N = 2000;  s = 2.0;
N = 5001;  s = 3.0;


[tx, ty, Nb] = rbfCenters.circleUniformCenters(N,1);

nh = Nb/2;
xc(1:nh) = tx(N-Nb+1:N-nh);        % half of boundary centers
xc(nh+1:N-nh) = tx(1:N-Nb);        % interior centers
xc(N-nh+1:N) = tx(N-nh+1:N);   % other half of boundary centers

yc(1:nh) = ty(N-Nb+1:N-nh);        % half of boundary centers
yc(nh+1:N-nh) = ty(1:N-Nb);        % interior centers
yc(N-nh+1:N) = ty(N-nh+1:N);   % other half of boundary centers

[xc,yc] = rbfCentro.centroCenters(xc,yc,2,false);
exact = 1 - xc + xc.*yc + 0.5*sin(pi*xc).*sin(pi*yc);
f = -pi^2*sin(pi*xc).*sin(pi*yc);

d = sqrt( xc.^2 + yc.^2 );
I = find( d<(1+100*eps) & d>(1-100*eps));  % locate boundary points
f(I) = exact(I);


if CENTRO
    
   tic
   [r, rx, ry] = phi.distanceMatrix2d(xc(1:N/2),yc(1:N/2),xc,yc);  % half-sized distance matrices
   B = phi.rbf(r,s);         % half-sized system matrix
   H = phi.L(r, s);          % Laplacian (half-sized)
   [kappaH, kappaL, kappaM] = rbfCentro.centroConditionNumber(H,mu)
    
    H(I,:) = B(I,:);  % Dirichlet BCs
    
    a = rbfCentro.solveCentro(H,f,mu,true); 
    [L,M] = rbfCentro.centroDecomposeMatrix(B,0);
    u = rbfCentro.centroMult(a,L,M,0);
    
    errorCentro = norm( u - exact, inf)
    toc
    
else
    
    tic
    [r, rx, ry] = phi.distanceMatrix2d(xc,yc);
    B = phi.rbf(r,s);
    H = phi.L(r, s);
    kappaH = cond(H)
    H(I,:) = B(I,:);
   
    a = H\f;
    u = B*a;
   
    errorStandard = norm( u - exact, inf)
    
    toc 
    
    rbfCentro.hasSymmetry(H);   % check the full sized matrix for symmetry

end



% scatter(xc,yc,'b.')
% hold on
% scatter(xc(I),yc(I),'ro')


% ---------- test for centrosymmetry with all boundary centers last -------
% -> the matrix is not centrosymmetric

% [xc, yc, Nb] = rbfCenters.circleUniformCenters(N,1);
% d = sqrt( xc.^2 + yc.^2 );
% I = find( d<(1+100*eps) & d>(1-100*eps));
% 
% [r, rx, ry] = phi.distanceMatrix2d(xc,yc);
% B = phi.rbf(r,s);
% H = phi.L(r, s);
% H(I,:) = B(I,:);
% rbfCentro.hasSymmetry(H); 




