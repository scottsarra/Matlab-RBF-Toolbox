% complexCentroCenters.m
% 
% Contructs a centro center distribution on a complexly shaped domain
% using quasi-random Hammersley points which are placed denser near the
% boundary than in the interior.  Before the centers are extended
% centrosymmetrically, the domain needs to be rotated so that it is symmetric
% with respect to the x-axis.  A different rotation could be used to make
% the domain symmetric with repect to the origin.

f = @(t) 0.8 + 0.1*( sin(6*t) + sin(3*t) );    % domain boundary

small = 0.005;
N = 5000;                  % N potiential centers in boundary region
boundaryLayerSize = 0.2;   % width of region with denser centers

% determine the size of the rectangle needed to cover the domain

t =  linspace(0,2*pi,200);          
x = f(t).*cos(t);  y = f(t).*sin(t);
A = min(x) - small; B = max(x) + small;
C = min(y) - small; D = max(y) + small;

[xc, yc] = rbfCenters.Hammersley2d(N);
xc = (B - A)*xc + A;             % [0,1] --> [A,B]
yc = (D - C)*yc + C;             % [0,1] --> [C,D]

% ---------- boundary region centers -----------------------

th = atan2(yc,xc);  p = sqrt(xc.^2 + yc.^2);
ro = f(th);                   % outter boundary
ri = ro - boundaryLayerSize;  % inner border of boundary region

xn = zeros(N,1);  yn = zeros(N,1);   I = 1;
for i=1:N
    if  and( p(i) < ro(i), p(i)>ri(i) )
        xn(I) = xc(i); yn(I) = yc(i);  I = I + 1;
    end
end
xc = xn(1:find(xn,1,'last'));    % remove trailing zeros
yc = yn(1:find(yn,1,'last'));

% ------ interior centers----------------------------

N2 = 1750;     % N2 potiential centers in interior region
[xci, yci] = rbfCenters.Hammersley2d(N2);

A = A + boundaryLayerSize;  B = B - boundaryLayerSize;
C = C + boundaryLayerSize;  D = D - boundaryLayerSize;
xci = (B - A)*xci + A;             % [0,1] --> [A,B]
yci = (D - C)*yci + C;             % [0,1] --> [C,D]

th = atan2(yci,xci); p = sqrt(xci.^2 + yci.^2);
ri = f(th) - boundaryLayerSize;    % interior region boundary

xn = zeros(N2,1);  yn = zeros(N2,1);  I = 1;
for i=1:N2
    if  p(i)<(ri(i) - small) 
        xn(I) = xci(i); yn(I) = yci(i);  I = I + 1;
    end
end
 
xci = xn(1:find(xn,1,'last'));    % remove trailing zeros
yci = yn(1:find(yn,1,'last'));
x = [xc; xci];   y = [yc; yci];   % merge centers

% --------- rotate the domain clockwise 0.25 radians ----------------------
% --- so that it is symmetric wrt the x-axis ------------------------------

t = atan2(y,x) - 0.25;  r = sqrt(x.^2 + y.^2);
x = r.*cos(t); y = r.*sin(t);

% --find centers in upper half of the domain (x-axis symmetry) ----------

I = find(y>(0 + 1e-3 )); x = x(I);  y = y(I);

% ------ extend "centrosymmetrically" to the other half -----------------

x = [x; flipud(x)];   y = [y; flipud(-y)];
[r, rx, ry] = rbfx.distanceMatrix2d(x,y);

disp(' ')        
% the signed distance matrix rx is NOT skew-centro => only even order DMs 
% wrt x will be centrosymmetric, odd order DMs wrt x will NOT be skew-centro
fprintf('rx: '); rbfCentro.hasSymmetry(rx);
% the signed distance matrix ry is skew-centro
fprintf('ry: '); rbfCentro.hasSymmetry(ry);
% the distance matrix r is centrosymmetric
fprintf('r: '); rbfCentro.hasSymmetry(r);    
 
% ---------------- rotate back to the original position -----------
% ------------ after any centro calculations are made -------------

 t = atan2(y,x) + 0.25;    r = sqrt(x.^2 + y.^2);
 x = r.*cos(t);  y = r.*sin(t);
             
% -----------------------------------------------------------------
scatter(x,y,'b.')
