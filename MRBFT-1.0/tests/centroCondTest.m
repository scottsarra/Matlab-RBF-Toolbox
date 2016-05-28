% centroCondTest
%
% Verifies the centrosymmetric condition number algorithm against the standard
% algorithm.  The two algorithms agree until the matrix becomes very 
% ill-conditioned.  As expected there is a slight variation when
% cond(B) > O(10^16)

warning off

phi = gax();
N = 44;
mu = 2e-16;  % MDI regularization parameter


%xc = linspace(-1,1,N)';
xc = -cos((0:N-1)*pi/(N-1))';    % centers

r = phi.distanceMatrix1d(xc(1:N/2),xc);     % left half of system matrix
rf = phi.distanceMatrix1d(xc,xc);           % full system matrix

sv = 10:-0.25:0.25;
Ns = length(sv);

cb = zeros(Ns,1); cbe = zeros(Ns,1);  cf = zeros(Ns,1);
for i=1:Ns
  s = sv(i);            
  B = phi.rbf(r,s);   % half-sized system matrix
  %  cbe(i) = phi.centroConditionNumberEig(B,mu); % ill-conditioning leads to complex condition number
  [kappaB, kappaL, kappaM] = rbfCentro.centroConditionNumber(B,mu); % the SVD version is more stable
  cb(i) = kappaB; 
  B = phi.rbf(rf,s);  % full system matrix
  cf(i) = cond(B + mu*eye(N));
end


semilogy(sv,cb,'b',sv,cf,'g--')
legend('cs kappa(B)','std kappa(B)')
warning on
