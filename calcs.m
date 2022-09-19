function P = calcs(price, theta,inp)

nj = inp.nj;
S = inp.S;
R = inp.R;
indivd = inp.indivd;
n = inp.n;

nss = inp.nss;
dkss = inp.dkss;
ntf = inp.ntf;
dytb = inp.dytb;
X = inp.prodcharac0;

% Chainging price
X(:,2) = price;

nx = size(X,2);
J = inp.J;

theta2 = theta(nx+1:end,1);
theta1 = theta(1:nx,1);

pincome = indivd(:,2);
age = indivd(:,4);
ij = inp.ij;


    
delta = X*theta1; % arbitrarily chosen starting values for delta

mu = 0; % util part related to theta2


% mu: j-by-R & nu:j-by-R
% xrc: price 
mu = mu + theta2(1)*((pincome)').*repmat(price,1,R) +theta2(2)*(age.').*repmat(nss,1,R);
%mu = mu + theta2(1)*((age)').*repmat(ntf,1,R) +theta2(2)*(age.').*repmat(ntf,1,R)+theta2(3)*(school.').*repmat(ntf,1,R)+theta2(4)*(plincome.').*repmat(ntf,1,R);
% theta2(end-1)*ppi'.*repmat(price,1,R),theta2(2)*(dhnet.').*repmat(ntf,1,R),theta2(4)*(female.').*repmat(nss,1,R)


U = repmat(delta,1,R) + mu;
eU = exp(U); %jm-by-R

denom0 = 1 +sum(eU,1);
denom = repmat(denom0,nj,1);

P = mean(eU./denom,2);

end

    




