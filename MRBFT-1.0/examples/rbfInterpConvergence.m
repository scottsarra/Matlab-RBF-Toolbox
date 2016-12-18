% % rbfInterpConvergence
%
% Convergence rate of a RBF interpolant with a fixed shape parameter and 
% increasing N.  The convergence is geometric (also called spectral or 
% exponential) as long as the floating point system can handle the condition 
% number of the system matrix.

warning off, clear, home, close all

phi = iqx();
mp.Digits(34);
mu = 0;        % MDI regularization parameter (no regularization)
safe = true;   % backslash rather than forcing Cholesky
s = mp('2.0');  pi = mp('pi');

M = 200;
x = linspace(mp(-1),mp(1),M)';             % evaluation points
fe = exp(sin(pi*x));

Nv = mp(5:10:110);
Ns = length(Nv);
er = mp( zeros(Ns,1) );

for k = 1:Ns
    
    N = Nv(k);
    
    xc = -cos(pi*mp(0:N-1)/(N-1))';   % boundary clustered centers
    r = rbfx.distanceMatrix1d(xc);
    f = exp(sin(pi*xc));
    
    B = phi.rbf(r,s);
    a = rbfx.solve(B,f,mu,safe);    % expansion coefficients

    re = rbfx.distanceMatrix1d(xc,x);
    H = phi.rbf(re,s);
    
    fa = H*a;
    er(k) = norm(fa - fe, inf);
end


semilogy(Nv,er,'b*')
xlabel('N'), ylabel('|error|')
warning on


figure()
loglog(er(1:end-2),er(2:end-1), 'g')

% rho approximately 1 implies geometric (spectral, or exponential) convergence
rho = (log10(er(end-1)) - log10(er(end-2)))/(log10(er(end-2)) - log10(er(end-3)))
