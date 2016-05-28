% mdiExample
%
% 1d interpolation problem using extended precision and regulatization by
% the method of diagonal increments.

warning off

phi = iqx();
safe = false;  % use Cholesky factorization
S = 1.45:-0.025:0.05;   % shape parameters
Sn = length(S);
M = 175;
mp.Digits(34);  N = mp('55');  pi = mp('pi');  
kappa = mp(zeros(Sn,1));  er = mp(zeros(Sn,1));  
kappa2 = mp(zeros(Sn,1));  er2 = mp(zeros(Sn,1));
I = mp( eye(N) );
mu = 10*mp.eps;      % MDI regularization parameter

gamma = 0.99;
xc = (asin(-gamma*cos(pi*(0:N-1)/(N-1)))/asin(gamma))'; % boundary clustered centers
r = rbfx.distanceMatrix1d(xc);
x = linspace(-1,1,M)';             % evaluation points
x = mp(x);
re = rbfx.distanceMatrix1d(xc,x);

f = exp(sin(pi*xc));
fe = exp(sin(pi*x));

for k = 1:Sn
    s = S(k);
    B = phi.rbf(r,s);
    kappa(k) = cond(B + mu*I);
    kappa2(k) = cond(B);
    a = rbfx.solve(B,f,mu,false);   % MDI and Cholesky factoization
    H = phi.rbf(re,s);
    fa = H*a;
    er(k) = norm(fa - fe, inf);
    
    a2 = rbfx.solve(B,f,0,true);   % no rugularization; use backslash
                                   % as Chol may fail without MDI
    fa2 = H*a2;
    er2(k) = norm(fa2 - fe, inf);
    
end

semilogy(S,kappa2,'b',S,kappa,'g--')
xlabel('shape parameter'), ylabel('\kappa(B)')

figure()

semilogy(S,er2,'b',S,er,'g--')
xlabel('shape parameter'), ylabel('|error|')


warning on