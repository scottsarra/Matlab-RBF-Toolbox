% mdiRegularization.m
%
% Interpolates the Franke function on scattered centers located in a domain
% that is one-fourth of a circle.  Compares the accuracy and condition number
% of the system matrix over a range of shape parameter with and without 
% regularization by the method of diagonal increments.

warning off, clear, home, close all
mu = 2e-14;

K = 1.5*sqrt(2);
% open text files located in the /examples folder
XC = dlmread('xc.txt',' ');  xc = XC(:,1)/K;  yc = XC(:,2)/K; % centers
 X = dlmread('x.txt',' ');    x = X(:,1)/K;    y = X(:,2)/K;  % evaluation points
     
N = length(xc); M = length(x);
     
fn = F2d();          %  Franke function
f = fn.F(xc,yc);
fe = fn.F(x,y);

phi = iqx();         % IQ RBF

[r, rx, ry] = rbfx.distanceMatrix2d(xc,yc);
[re, rx, ry] = rbfx.distanceMatrix2d(xc,yc,x,y);
    
S = 6:-0.1:0.2;     Sn = length(S);
kappa = zeros(Sn,1);  er = zeros(Sn,1);  
kappa2 = zeros(Sn,1);  er2 = zeros(Sn,1);
I = eye(N);

for k = 1:Sn
    s = S(k);
    B = phi.rbf(r,s);  % system matrix

    kappa(k) = cond(B + mu*I);
    kappa2(k) = cond(B);
    a = rbfx.solve(B,f,mu,false);   % regularized system solver
   
    H = phi.rbf(re,s);
    fa = H*a;
    er(k) = norm(fa - fe, inf);
    
    a2 = rbfx.solve(B,f,0,true);   % no regularization
    fa2 = H*a2;
    er2(k) = norm(fa2 - fe, inf);
  
  end
 
warning on

semilogy(S,kappa2,'b--',S,kappa,'g')
xlabel('shape parameter'), ylabel('\kappa(B)')

figure()

semilogy(S,er2,'b--',S,er,'g')
xlabel('shape parameter'), ylabel('|error|')