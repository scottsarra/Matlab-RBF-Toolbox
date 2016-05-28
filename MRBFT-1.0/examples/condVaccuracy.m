% condVaccuracy.m

warning off
QUADRUPLE = false;

phi = iqx();
mu = 0;           % MDI regularization parameter
safe = true;      % backslash rather than forcing Cholesky
S = 3:-0.1:0.1;   % shape parameters
Sn = length(S);
M = 175;

if QUADRUPLE   
   mp.Digits(34);  N = mp('55');  pi = mp('pi');  kappa = mp(zeros(Sn,1));  er = mp(zeros(Sn,1));  % quadruple  
else
    N = 55; kappa = zeros(Sn,1);  er = zeros(Sn,1);  % double 
end


gamma = 0.99;
xc = (asin(-gamma*cos(pi*(0:N-1)/(N-1)))/asin(gamma))';   % boundary clustered centers
r = rbfx.distanceMatrix1d(xc);
x = linspace(-1,1,M)';                                    % evaluation points
x = mp(x);
re = rbfx.distanceMatrix1d(xc,x);

f = exp(sin(pi*xc));
fe = exp(sin(pi*x));

for k = 1:Sn
    s = S(k);
    B = phi.rbf(r,s);
    kappa(k) = cond(B);
    a = rbfx.solve(B,f,mu,safe);
    H = phi.rbf(re,s);
    fa = H*a;
    er(k) = norm(fa - fe, inf);
end

%semilogy(S,kappa,'b')    % plot condition number versus shape
%xlabel('shape parameter'), ylabel('\kappa(B)')
% 
% figure()

semilogy(S,er,'b')       % plot error versus shape parameter
xlabel('shape parameter'), ylabel('|error|')

warning on