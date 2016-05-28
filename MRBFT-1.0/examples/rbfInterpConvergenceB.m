% rbfInterpConvergenceB
%
% Similiar to rbfInterpConvergence.m except that the number of centers N is
% fixed and the shape parameter is decreasing. The convergence exponential 
% as long as the floating point system can handle the condition 
% number of the system matrix.

warning off, clear, home, close all


phi = iqx();
mp.Digits(34);
mu = 0;         % MDI regularization parameter (no regularization)
safe = true;    % backslash rather than forcing Cholesky
s = mp('2.0');  pi = mp('pi');

N = mp('90');
xc = -cos(pi*(0:N-1)/(N-1))'; % boundary clustered centers
r = rbfx.distanceMatrix1d(xc);
f = exp(sin(pi*xc));

M = 200;
x = mp( linspace(-1,1,M)' );             % evaluation points
fe = exp(sin(pi*x));

re = rbfx.distanceMatrix1d(xc,x);
 
Sv = mp('10'):mp('-0.25'):mp('2.5');
Ns = length(Sv);
er = mp( zeros(Ns,1) );

for k = 1:Ns
    s = Sv(k);
    B = phi.rbf(r,s);
    a = rbfx.solve(B,f,mu,safe);
    H = phi.rbf(re,s);
    fa = H*a;
    er(k) = norm(fa - fe, inf);
end


semilogy(Sv,er,'b*')
xlabel('shape parameter'), ylabel('|error|')
warning on


figure()
loglog(er(1:end-2),er(2:end-1), 'g')

% rho approximately 1 implies exponential convergence
rho = (log10(er(end-1)) - log10(er(end-2)))/(log10(er(end-2)) - log10(er(end-3)))