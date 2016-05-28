% diffusionReactionCentroDriver.m

% Solves the diffusion-reaction PDE
%  u_t  = visc*(u_xx + u_yy) + gamma*u^2(1-u)
% on a circular domain with Dirichlet boundary conditions prescribed from
% the exact solution.  The PDE is discretized in space with the IQ RBF method
% and is advanced in time with a 4th order Runge-Kutta method.
%
% The problem is solved two ways: 1) standard algorithms, 2) centrosymmetric 
% algorithms.  With N = 5000 the accuracy of the 2 approaches is the same
% but the centrosymmetric approach is approximately 8 times faster and 
% requires half the storage.

clear, home, format compact

visc = 0.5;    % viscosity coefficient
dt = 0.001;    % time step size
finalT = 5;    % advance in time from t = 0 to t = finalT


CENTRO = false

%N = 2000;  shape = 0.4;
N = 5001;  shape = 0.4;

[tx, ty, Nb] = rbfCenters.circleUniformCenters(N,1);

nh = Nb/2;
xc(1:nh) = tx(N-Nb+1:N-nh);        % half of boundary centers
xc(nh+1:N-nh) = tx(1:N-Nb);        % interior centers
xc(N-nh+1:N) = tx(N-nh+1:N);   % other half of boundary centers

yc(1:nh) = ty(N-Nb+1:N-nh);        % half of boundary centers
yc(nh+1:N-nh) = ty(1:N-Nb);        % interior centers
yc(N-nh+1:N) = ty(N-nh+1:N);   % other half of boundary centers

[xc,yc] = rbfCentro.centroCenters(xc,yc,2,false);

d = sqrt( xc.^2 + yc.^2 );                % locate boundary centers
I = find( d<(1+100*eps) & d>(1-100*eps));


xc = 10*xc + 5;   % domain of a radius 10 circle centered at (5,5)
yc = 10*yc + 5;

sol = diffusionReactionCentro(xc,yc,I,visc,dt,finalT,CENTRO,shape);

td = delaunay(xc,yc);
trisurf(td,xc,yc,sol);