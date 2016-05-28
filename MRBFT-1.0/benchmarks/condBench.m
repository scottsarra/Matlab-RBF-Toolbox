% condBench.m
%
% Compares the execution times of centrosymmetric versus standard algorithms for
%   calculating the 2-norm condition number of a centrosymmetric matrix

clear, home, format compact

phi = iqx();
s = 230;
mu = 0;
its = 5;

Nv = 600:250:2100;

Ns = length(Nv);
cTime = zeros(Ns,1);
ncTime = zeros(Ns,1);


tic
for k=1:Ns
    N = Nv(k);
    x = linspace(-100,100,N)';
    f = rand(N,1);
    tic
    for j=1:its
        r = phi.distanceMatrix1d(x(1:N/2),x);  % construct half-sized distance matrix
        B = phi.rbf(r,s);                      % half-sized system matrix
        [kappaB, kappaL, kappaM] = rbfCentro.centroConditionNumber(B,mu);
    end
    cTime(k) = toc;
end
    
    
for k=1:Ns
    N = Nv(k);
    x = linspace(-1,1,N)';
    f = rand(N,1);
    tic
    for j=1:its
        r = phi.distanceMatrix1d(x);
        B = phi.rbf(r,s);
        kappa = cond(B);
    end
    ncTime(k) = toc;
end


plot(Nv,cTime./ncTime,'b*',Nv,ones(Ns,1),'g--')
xlabel('N'), ylabel('execution time ratio')