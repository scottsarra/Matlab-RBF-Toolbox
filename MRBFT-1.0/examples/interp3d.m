% interp3d.m
%
% Gaussian RBF interpolation on the surface of a sphere.
%

warning off
tic

phi = gax();
mu = 1.5e-13;
safe = false;
N = 6000;
M = 7000;

% http://www.mathworks.com/matlabcentral/fileexchange/6977-pointonsphere/content/pointonsphere.m
P = pointonsphere(N);
xc = P(:,1); yc = P(:,2); zc = P(:,3);

P = pointonsphere(M);
x = P(:,1); y = P(:,2); z = P(:,3);

f = 0.1*( 9*xc.^3 - 2*xc.^2.*yc + 3*xc.*yc.^2 - 4*yc.^3 + 2*zc.^3 - xc.*yc.*zc );
fe = 0.1*( 9*x.^3 - 2*x.^2.*y + 3*x.*y.^2 - 4*y.^3 + 2*z.^3 - x.*y.*z );

[r,rx,ry,rz] = phi.distanceMatrix3d(xc,yc,zc);
[re,rx,ry,rz] = phi.distanceMatrix3d(xc,yc,zc,x,y,z);


sv = 8.1:-1:0.1;
Ns = length(sv);

er = zeros(Ns,1);
for i=1:Ns
    s = sv(i)      
    B = phi.rbf(r,s);
    a = phi.solve(B,f,mu,safe);
    H = phi.rbf(re,s);
    fa = H*a;
    er(i) = norm(fa - fe, inf);
end

toc
semilogy(sv,er,'b--')
warning
