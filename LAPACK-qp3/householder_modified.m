function [tau,v,beta] = householder_modified(x)
if nargin < 1,
  error('Reflector requires at least one argument');
end
if nargin < 2,
   j = 1; 
end
if nargin < 3,
   % Ideally this should be safeMin/epsilon
   safeInv = 2e-292;
end

alpha=x(j);
nrmExl=norm(x([1:j-1 j+1:end]));

v=x;
if nrmExl == 0 && imag(alpha) == 0,
  beta=-x(j);
  v(j)=1;
  tau=2;
  return;
end
if real(alpha) <= 0,
  beta = norm([alpha;nrmExl]);
else 
  beta = -norm([alpha;nrmExl]);
end

% Rescale if necessary
count=0;
if abs(beta) < safeInv,
  invOfSafeInv = 1/safeInv;
  while abs(beta) < safeInv,
    count=count+1;
    v([1:j-1 j+1:end])=x([1:j-1 j+1:end])*invOfSafeInv;
    alpha = alpha*invOfSafeInv;
    beta = beta*invOfSafeInv;
  end  
  nrmExl = norm(v([1:j-1 j+1:end]));
  if real(alpha) <= 0,
    beta = norm([alpha;nrmExl]);
  else
    beta = -norm([alpha;nrmExl]);
  end
end

tau = (beta-conj(alpha))/beta;
v(j)=1;
v([1:j-1 j+1:end])=v([1:j-1 j+1:end])/(alpha-beta);

% Undo the scaling
for k=1:count,
  beta = beta*safeInv;
end