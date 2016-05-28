% systemSolveBench.m
%
% Compares the evaluation times of centrosymmetric versus standard algorithms
% for the solution of a centrosymmetric linear system

clear, home, format compact

phi = iqx();
s = 330;
rho = 0;
safe = false;  % use Cholesky factorization
mu = 5e-15;    % MDI regularization parameter
its = 5;

Nv = 350:250:4100;
Ns = length(Nv);
cTime = zeros(Ns,1);
ncTime = zeros(Ns,1);

tic
for k=1:Ns     % centrosymmetric 
    N = Nv(k);
    x = linspace(-100,100,N)';
    f = rand(N,1);
    tic
    for j=1:its
        r = phi.distanceMatrix1d(x(1:N/2),x);    % half sized distance matrix
        B = phi.rbf(r,s);                        % half sized system matrix
        a = rbfCentro.solveCentro(B,f,mu,safe);
    end
    cTime(k) = toc;
end
       
for k=1:Ns    % standard
    N = Nv(k);
    x = linspace(-1,1,N)';
    f = rand(N,1);
    tic
    for j=1:its
        r = phi.distanceMatrix1d(x);
        B = phi.rbf(r,s);
        a = phi.solve(B,f,mu,safe);
    end
    ncTime(k) = toc;
end

plot(Nv,cTime./ncTime,'b*')
xlabel('N'), ylabel('execution time ratio')