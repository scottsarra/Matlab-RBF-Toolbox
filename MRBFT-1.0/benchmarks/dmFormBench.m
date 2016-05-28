% dmFormBench.m
%
% compares the execution times of centro versus standard algorithms
%   for constructing a centro derivative matrix

clear, home, format compact

phi = iqx();
s = 230;
rho = 1;
safe = true;
mu = 0;
its = 10;

Nv = 350:250:2100;
Ns = length(Nv);
cTime = zeros(Ns,1);
ncTime = zeros(Ns,1);


tic
for k=1:Ns
    N = Nv(k);
    x = linspace(-1,1,N)';
    tic
    for j=1:its
        r = phi.distanceMatrix1d(x(1:N/2),x);  % half-sized distance matrix
        B = phi.rbf(r,s);
        H = phi.D1(r,s,r);
        D = rbfCentro.centroDM(B,H,N,rho,mu,safe);
    end
    cTime(k) = toc;
end
    
    
for k=1:Ns
    N = Nv(k);
    x = linspace(-1,1,N)';
    tic
    for j=1:its
        r =  phi.distanceMatrix1d(x);  % full-sized distance matrix
        B = phi.rbf(r,s);
        H = phi.D1(r,s,r);
        D = phi.dm(B,H,mu,safe);
    end
    ncTime(k) = toc;
end


plot(Nv,cTime./ncTime,'b*')
xlabel('N'), ylabel('execution time ratio')

