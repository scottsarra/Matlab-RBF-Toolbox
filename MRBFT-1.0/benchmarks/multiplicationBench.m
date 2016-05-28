% multiplicationBench.m
%
% compares the execution time of centro versus the standard algorithm for
%   matrix-vector multiplication

clear, home, format compact

phi = iqx();
s = 230;
rho = 1;
safe = true;  % use backslash operator
mu = 0;       % no regularization
its = 25;

Nv = 600:100:3000;
Ns = length(Nv);
cTime = zeros(Ns,1);
ncTime = zeros(Ns,1);


for k=1:Ns
    N = Nv(k);
    x = linspace(-1,1,N)';
    u = rand(N,1);
    r = phi.distanceMatrix1d(x(1:N/2),x);
    B = phi.rbf(r,s);
    H = phi.D1(r,s,r);
    D = rbfCentro.centroDM(B,H,N,rho,mu,safe);
    [L,M] = rbfCentro.centroDecomposeMatrix(D,rho);
    tic
    for j=1:its   
        uac = rbfCentro.centroMult(u,L,M,rho);
    end
    cTime(k) = toc;
end
    
    
for k=1:Ns
    N = Nv(k);
    u = rand(N,1);
    x = linspace(-1,1,N)';
    r =  phi.distanceMatrix1d(x);
    B = phi.rbf(r,s);
    H = phi.D1(r,s,r);
    D = phi.dm(B,H,mu,safe);
    tic
    for j=1:its
        ua = D*u;
    end
    ncTime(k) = toc;
end

plot(Nv,cTime./ncTime,'b*')
xlabel('N'), ylabel('execution time ratio')

