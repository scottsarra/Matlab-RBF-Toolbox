% isCentroTest.m
%
% Depending on how the centers were extended to be symmetric, RBF
% differentiation matrices will have a (skew) centrosymmetric structure.
% The following reference can be consulted for details:
% "Radial Basis Function Methods - the case of symmetric domains."  
% Under review, Numerical Methods for Partial Differential Equations, 2016.

phi = iqx();

s = 10;  % shape parameter

%  symType     0 - y axis
%              1 - x axis  
%              2 - origin  -> all order derivative have correct symmetry

symType = 2;

%[xc,yc] = rbfCentro.centroCircle(500,true,0,1,false);   % uses origin sym


[x, y] = rbfCenters.circleCenters(500,true,0,1,false);  
[xc,yc] = rbfCentro.centroCenters(x,y,symType,true);             


[r, rx, ry] = rbfx.distanceMatrix2d(xc,yc);

H1 = phi.D1(r,s,rx);
H2 = phi.D2(r,s,rx);
H3 = phi.D3(r,s,rx);
H4 = phi.D4(r,s,rx);
G = phi.G(r, s, rx, ry);
L = phi.L(r, s);
B = phi.B(r, s, rx, ry) ;
H12 = phi.D12(r, s, rx, ry);
H22 = phi.D22(r, s, rx, ry);

disp(' ')
fprintf('D1: '); rbfCentro.hasSymmetry(H1);
fprintf('D2: '); rbfCentro.hasSymmetry(H2);
fprintf('D3: '); rbfCentro.hasSymmetry(H3);
fprintf('D4: '); rbfCentro.hasSymmetry(H4);

fprintf('gradient: '); rbfCentro.hasSymmetry(G);
fprintf('Laplacian: '); rbfCentro.hasSymmetry(L);
fprintf('Biharmonic: '); rbfCentro.hasSymmetry(B);

fprintf('mixed partial 12: '); rbfCentro.hasSymmetry(H12);
fprintf('mixed partial 22: '); rbfCentro.hasSymmetry(H2);