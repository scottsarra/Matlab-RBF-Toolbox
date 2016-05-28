% rbfDerivativeTest
%
% Tests all derivative approximation methods of the iqx and gax classes
%   using both double and quadruple precision


% ---------- output of the four tests -------------------  


% OUTPUT OF TEST 1
%dx error = 6.83e-09
%dxx error = 1.75e-07
%dxxx error = 8.82e-06
%dxxxx error = 6.47e-04
%dy error = 7.34e-09
%dyy error = 1.24e-07
%dyyy error = 8.73e-06
%dyyyy error = 8.38e-04
%gradient error = 1.05e-08
%Laplacian error = 1.68e-07
%biharmonic error = 1.36e-03
%dx1y2 error = 4.60e-06
%dx2y1 error = 4.43e-06
%dx2y2 error = 2.94e-04

% ---------------------------

% OUTPUT OF TEST 2
%dx error = 9.44e-12
%dxx error = 3.03e-10
%dxxx error = 2.12e-08
%dxxxx error = 5.31e-06
%dy error = 1.06e-11
%dyy error = 3.27e-10
%dyyy error = 6.19e-08
%dyyyy error = 1.27e-05
%gradient error = 1.86e-11
%Laplacian error = 2.74e-10
%biharmonic error = 1.95e-05
%dx1y2 error = 3.00e-08
%dx2y1 error = 2.13e-08
%dx2y2 error = 4.23e-06

% ---------------------------- 

% OUTPUT OF TEST 3
%dx error = 3.40e-05
%dxx error = 5.92e-03
%dxxx error = 7.53e-01
%dxxxx error = 7.11e+01
%dy error = 4.30e-05
%dyy error = 6.07e-03
%dyyy error = 7.38e-01
%dyyyy error = 6.68e+01
%gradient error = 5.74e-05
%Laplacian error = 6.86e-03
%biharmonic error = 1.14e+02
%dx1y2 error = 3.99e-01
%dx2y1 error = 3.00e-01
%dx2y2 error = 2.55e+01

% ---------------------------

% OUTPUT OF TEST 4
%dx error = 1.04e-04
%dxx error = 2.79e-03
%dxxx error = 3.40e-01
%dxxxx error = 4.01e+01
%dy error = 9.49e-05
%dyy error = 3.11e-03
%dyyy error = 3.77e-01
%dyyyy error = 2.98e+01
%gradient error = 1.40e-04
%Laplacian error = 3.25e-03
%biharmonic error = 6.00e+01
%dx1y2 error = 1.32e-01
%dx2y1 error = 1.91e-01
%dx2y2 error = 1.34e+01

% ----------------------------------------------------------------------

clear, home, format compact

TESTNUMBER = 4;

if TESTNUMBER == 1
    % test 1, IQ with quadruple precision
    phi = iqx();
    mp.Digits(34);  s = mp('1.2'); N = mp('2000'); 
elseif TESTNUMBER == 2
    % test 2, GA with quadruple precision
    phi = gax();
    mp.Digits(34);  s = mp('3.5'); N = mp('2000');    
elseif TESTNUMBER == 3
    % test 3, GA with double precision
    phi = gax();
    s = 4.5; N = 2000; 
else
    % test 4, IQ with double precision
    phi = iqx();
    s = 2.35; N = 2000; 
end

nCh = inf;                  % norm choice
G = F2c;

[x, y] = rbfCenters.circleCenters(N,true,1,1,false);
f = G.F(x,y);

[r, rx, ry] = phi.distanceMatrix2d(x,y);
B = phi.rbf(r,s);

mu = 0;
safe = true;  % use mldivide rather than Cholesky directly
a = rbfx.solve(B,f,mu,safe);


H = phi.D1(r,s,rx);
fx = H*a;
fprintf('dx error = %4.2e\n',norm(fx - G.x1(x,y), nCh));

H = phi.D2(r,s,rx);
fxx = H*a;
fprintf('dxx error = %4.2e\n',norm(fxx - G.x2(x,y), nCh));

H = phi.D3(r,s,rx);
fxxx = H*a;
fprintf('dxxx error = %4.2e\n',norm(fxxx - G.x3(x,y), nCh));

H = phi.D4(r,s,rx);
fxxxx = H*a;
fprintf('dxxxx error = %4.2e\n',norm(fxxxx - G.x4(x,y), nCh));


H = phi.D1(r,s,ry);
fy = H*a;
fprintf('dy error = %4.2e\n',norm(fy - G.y1(x,y), nCh));

H = phi.D2(r,s,ry);
fyy = H*a;
fprintf('dyy error = %4.2e\n',norm(fyy - G.y2(x,y), nCh));

H = phi.D3(r,s,ry);
fyyy = H*a;
fprintf('dyyy error = %4.2e\n',norm(fyyy - G.y3(x,y), nCh));

H = phi.D4(r,s,ry);
fyyyy = H*a;
fprintf('dyyyy error = %4.2e\n',norm(fyyyy - G.y4(x,y), nCh));


H = phi.G(r,s,rx, ry);
fG = H*a;
fprintf('gradient error = %4.2e\n',norm(fG - G.G(x,y), nCh));

H = phi.L(r,s);
fL = H*a;
fprintf('Laplacian error = %4.2e\n',norm(fL- G.L(x,y), nCh));

H = phi.B(r,s,rx, ry);
fB = H*a;
fprintf('biharmonic error = %4.2e\n',norm(fB - G.B(x,y), nCh));

H = phi.D12(r,s,rx,ry);
f12 = H*a;
fprintf('dx1y2 error = %4.2e\n',norm(f12- G.p12(x,y), nCh));

H = phi.D12(r,s,ry,rx);
f21 = H*a;
fprintf('dx2y1 error = %4.2e\n',norm(f21- G.p21(x,y), nCh));

H = phi.D22(r,s,ry,rx);
f22 = H*a;
fprintf('dx2y2 error = %4.2e\n',norm(f22- G.p22(x,y), nCh));