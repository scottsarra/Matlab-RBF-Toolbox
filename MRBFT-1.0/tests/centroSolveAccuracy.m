% centroSolveAccuracy.m
%
% Compares the accuracy of centrosymmetric versus standard algorithms for solving
%  a centrosymmetic linear system.  The linear system is the system for the RBF
%  expansion coefficients over a range of the shape parameter. The centrosymmetric 
%  algorithm is slightly more accurate at most shape parameters and several 
%  decimal places more accurate for several shape parameters.

warning off
phi = gax();
N = 44;
mu = 1e-15;    % regularization parameter
safe = false;  % use Cholesky factorization

%xc = linspace(-1,1,N)';
xc = -cos((0:N-1)*pi/(N-1))';    % boundary clustered CGL centers

r = phi.distanceMatrix1d(xc(1:N/2),xc);     % left half of system matrix
rf = phi.distanceMatrix1d(xc,xc);           % full system matrix

sv = 12:-0.25:3.0;
Ns = length(sv);
o = ones(N,1);          % exact solution of the linear system

er = zeros(Ns,1);  er2 = zeros(Ns,1);  
for i=1:Ns
  s = sv(i);   
  
  B = phi.rbf(r,s);                         % half-sized system matrix
  [L,M] = rbfCentro.centroDecomposeMatrix(B,0); 
  f = rbfCentro.centroMult(o,L,M,0);              % f = B*o, right side so that o is the exact solution
  a = rbfCentro.solveCentro(B,f,mu,safe);
  er(i) = norm(a - o, inf);                 % error from centrosymmetric solver
  
  B2 = phi.rbf(rf,s);                       % full-sized system matrix
  f = B2*o;
  a2 = phi.solve(B2,f,mu,safe);
  er2(i) = norm(a2 - o, inf);               % error from solving the full system
  
end

semilogy(sv,er,'g*',sv,er2,'b*')
legend('centrosymmetric','standard')
xlabel('shape parameter'), ylabel('|error|')
warning







































% H NOT CENTRO
%Nv = 8:2:16;
%mu = 5e-14;
%safe = false;
%
%Ns = length(Nv);
%er = zeros(Ns,1);  erc = zeros(Ns,1); 
%
%for k = 1:Ns
%    N = Nv(k);
%    ae = ones(N,1);
%    H = hilb(N);
%    f = H*ae;
%    ac = rbfx.solveCentro(H,f,mu,safe);
%    erc(i) = norm(ac - ae, inf); 
%    a = rbfx.solve(H,f,mu,safe);
%    er(i) = norm(a - ae, inf); 
%end
%
%semilogy(Nv,er,'b',Nv,erc,'g')

%warning off
%tic
%
%phi = gax();
%N = 44;
%M = 175;
%mu = 2e-15;
%safe = false;
%
%xc = linspace(-1,1,N)';
%%xc = -cos((0:N-1)*pi/(N-1))';    % centers
% x = linspace(-1,1,M)';
%
%%f = exp(sin(pi*xc));
%%fe = exp(sin(pi*x));
%
%func = F1a();
%f = func.F(xc);
%fe = func.F(x);
%
%r = phi.distanceMatrix1d(xc(1:N/2),xc);     % left half of system matrix
%re = phi.distanceMatrix1d(xc,x);
%
%rf = phi.distanceMatrix1d(xc,xc);           %  full system matrix
%
%
%sv = 2.5:-0.01:2.0;
%Ns = length(sv);
%
%er = zeros(Ns,1);  er2 = zeros(Ns,1);  erL = zeros(Ns,1);  erM = zeros(Ns,1);
%for i=1:Ns
%  s = sv(i);            
%  B = phi.rbf(r,s);
%  a = phi.solveCentro(B,f,mu,safe);
%  H = phi.rbf(re,s);
%  fa = H*a;
%  er(i) = norm(fa - fe, inf);                 % interpolation error
%  %[kappaB, kappaL, kappaM, ews] = phi.centroSystemMatrixCond(B,mu);       % system matrix condition number
%  %er(i) = kappaB; erL(i) = kappaL;  erM(i) = kappaM;
%  
%%  B = phi.fullCentroMatrix(B,N,0);
%  B2 = phi.rbf(rf,s);
%  a2 = phi.solve(B2,f,mu,safe);
%  fa2 = H*a2;
%  er2(i) = norm(fa2 - fe, inf); 
%  
%%  er2(i) = cond(B + mu*eye(N));
%end
%
%toc
%semilogy(sv,er,'g',sv,er2,'b')
%warning
%min(er)