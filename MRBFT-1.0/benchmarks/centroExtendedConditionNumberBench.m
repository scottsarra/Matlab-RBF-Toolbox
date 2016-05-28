% centroExtendedConditionNumberBench.m
%
% Compares the execution times of centrosymmetric and standard algorithms
% for the 2-norm condition number using both double and quadruple precision.

% sample output (old MCT on Linux)
% centroExtendedConditionNumberBench
% double: standard = 0.574
% double: centro = 0.149
% double: standard to centro = 3.850
%  
% extended: standard = 154.553
% extended: centro = 51.310
% extended: standard to centro = 3.012

clear, home, format compact

phi = iqx();

N = 2000;
s = 8;
mu = 0;
[xc,yc] = rbfCentro.centroCircle(N,true,0,1,false);
N = length(xc);

% ------------------- double ------------------------------------

tic
r = phi.distanceMatrix2d(xc(1:N/2),yc(1:N/2),xc,yc);  % construct half-sized distance matrix
B = phi.rbf(r,s);                      % half-sized system matrix
[kappaB, kappaL, kappaM] = rbfCentro.centroConditionNumber(B,mu);
cTime = toc;

    

tic
r = phi.distanceMatrix2d(xc,yc);
B = phi.rbf(r,s);
kappa = cond(B);
ncTime = toc;


fprintf('double: standard = %4.3f\n',ncTime);
fprintf('double: centro = %4.3f\n',cTime);
fprintf('double: standard to centro = %4.3f\n',ncTime/cTime);
disp(' ')


% -------------------- extended ---------------------------------------

mp.Digits(34); 

N = 2000;
s = mp('8');
mu = 0;
[xc,yc] = rbfCentro.centroCircle(N,true,0,1,false);
N = length(xc);
xc = mp(xc);  yc = mp(yc);

% ------------------------------------------------------------------------

tic
r = phi.distanceMatrix2d(xc(1:N/2),yc(1:N/2),xc,yc);  % construct half-sized distance matrix
B = phi.rbf(r,s);                      % half-sized system matrix
[kappaB, kappaL, kappaM] = rbfCentro.centroConditionNumber(B,mu);
cTime = toc;

    

tic
r = phi.distanceMatrix2d(xc,yc);
B = phi.rbf(r,s);
kappa = cond(B);
ncTime = toc;


fprintf('extended: standard = %4.3f\n',ncTime);
fprintf('extended: centro = %4.3f\n',cTime);
fprintf('extended: standard to centro = %4.3f\n',ncTime/cTime);
disp(' ')


