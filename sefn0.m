%This calculates variace of blp estimates
function [gmmvar] = sefn0(theta,inp)

nj = inp.nj;
Z = inp.Z;

% m: total number of moments == total number of instruments
m = size(Z,2);

% Gamma
step = 1e-5; % step-length for forward-backward finite differences
K = length(theta);


Gamma = zeros(m,K); 

Dgni = zeros(nj,m,2);
for k=1:K
    for t=1:2 % 1=backward step; 2=forward step;
        theta_eps = theta; % "theta + epsilon" (small perturbation)
        theta_eps(k) = theta_eps(k) + ((-1)^t)*step;
        Dgni(:,:,t) = gi(theta_eps,inp);
    end
    Dgn = sum(Dgni,1);
    Gamma(:,k) = (diff(Dgn,1,3)')/(2*step); % diff between second col. and first col.
end
Gamma = (1/nj).*Gamma;

% V
gni = gi(theta,inp);

gni_squared = zeros(m,m,nj);

for i =1:nj
    gni_squared(:,:,i) = (gni(i,:)')*gni(i,:);
end
V = ((1/nj)^2).*sum(gni_squared,3);

invv = inv(V);
A = (1/nj).*cov(gni);

W = inv(A);




m1 = Gamma'*W*Gamma;
m2 = Gamma'*invv*Gamma;
m3 = inv(m1);
m4 = inv(m2);

% gmmvar0 = inv(m1*(m2\m1));
% gmmvar1 = inv(Gamma'*W*Gamma)*(Gamma'*W*V*W*Gamma)*inv(Gamma'*W*Gamma) ;
gmmvar = inv(Gamma'*(V\Gamma));
end