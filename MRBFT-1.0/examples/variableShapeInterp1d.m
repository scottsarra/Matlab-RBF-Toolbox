% variableShapeInterp1d 

% Variable shape parameter versus constant shape. This is a typical example
% in which the two approaches have system matrices with approximately the  
% same condition number, but the variable shape approach is several decimal
% places more accurate.


clear, home, warning off, format compact

phi = gax();
N = 44;  M = 175;
safe = true;       

%xc = linspace(-1,1,N)';
xc = -cos((0:N-1)*pi/(N-1))';    % centers
 x = linspace(-1,1,M)';

% problem 1: a very smooth function
% f = exp(sin(pi*xc));  fe = exp(sin(pi*x)); fp = pi*cos(pi*xc).*exp(sin(pi*xc));

% problem 2: a constant function
 f = ones(N,1);   fe = ones(M,1);   fp = zeros(N,1);

r = phi.distanceMatrix1d(xc);
re = phi.distanceMatrix1d(xc,x);


% ------------- constant shape, interpolation -------------------

s = 1.0;       
B = phi.rbf(r,s);
kappa = cond(B);
a = phi.solve(B,f,0,safe); 
H = phi.rbf(re,s);
fa = H*a;
er = norm(fa - fe, inf);


% ------------- variable shape, interpolation -------------------

sMin = 0.5;
sMax = 1.5;
opt = 3;

[sn, sm] = phi.variableShape(sMin,sMax,opt,N,M); 

Bv = phi.rbf(r,sn);
kappaV = cond(Bv);
av = phi.solve(Bv,f,0,safe); 
Hv = phi.rbf(re,sm);
fav = Hv*av;
erV = norm(fav - fe, inf);


% ------- constant shape, derivative  ----------------------------


H = phi.D1(r,s,r);
fa = H*a;
erd = norm(fa - fp, inf);

% ------- variable shape, derivative  ----------------------------

Hv = phi.D1(r,sn,r);
fav = Hv*av;
erdv = norm(fav - fp, inf);

fprintf('variable shape system matrix condition number: %1.2e \n',kappaV)
fprintf('constant shape system matrix condition number: %1.2e \n\n',kappa)

fprintf('variable shape interpolation error: %1.2e \n',erV)
fprintf('constant shape interpolation error: %1.2e \n\n',er)


fprintf('variable shape derivative error: %1.2e \n',erdv)
fprintf('constant shape derivative error: %1.2e \n',erd)