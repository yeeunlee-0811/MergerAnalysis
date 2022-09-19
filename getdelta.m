% this function gives deltas and choice probabilities P
function [delta, P] = getdelta(theta2,inp)

nj = inp.nj;
S = inp.S;
R = inp.R;
indivd = inp.indivd;


logS = log(S);

price = inp.price;
nss = inp.nss;
dkss = inp.dkss;
ntf = inp.ntf;
dytb = inp.dytb;
ihhnet = inp.ihhnet;

pincome = indivd(:,2);
age = indivd(:,4);
school = inp.school;
female = inp.female;
dmar = inp.dmar;
djob = inp.djob;
dytub = inp.dytub;
payer = inp.payer;
pi = inp.pi;
nu = inp.nu;
netflix = inp.netflix;
koreanss = inp.koreanss;
youtube = inp.youtube;

nu1 = nu(:,1);
nu2 = nu(:,2);
nu3 = nu(:,3);

nuy = nu1.*youtube;
nun = nu2.*netflix;
nuk = nu3.*koreanss;



%%%% First without random coefficient
% xrc = prodcharac;

tol = 1e-5;
dst = 1;
delta = 0.5*ones(nj,1); % arbitrarily chosen starting values for delta

mu = 0; % util part related to theta2


% mu: j-by-R & nu:j-by-R
% xrc: price 
mu = mu + theta2(1)*(pincome.').*repmat(price,1,R)+theta2(2)*(age.').*repmat(nss,1,R)+theta2(3)*(ihhnet.').*repmat(ntf,1,R);
%theta2(1)*(ihhnet.').*repmat(ntf,1,R)
%mu = mu +theta2(1)*(nuy.').*repmat(dytb,1,R)+theta2(2)*(nun.').*repmat(ntf,1,R)+theta2(3)*(nuk.').*repmat(dkss,1,R);
%theta2(1)*(pincome.').*repmat(price,1,R)+theta2(2)*(age.').*repmat(nss,1,R)
%mu = mu + theta2(1)*(pincome.').*repmat(price,1,R) +theta2(2)*(pi.').*repmat(nss,1,R)+theta2(3)*(youtube.').*repmat(dytb,1,R)+theta2(3)*(netflix.').*repmat(ntf,1,R)+theta2(4)*(koreanss.').*repmat(dkss,1,R);
%mu = mu + theta2(1)*((age)').*repmat(ntf,1,R) +theta2(2)*(age.').*repmat(ntf,1,R)+theta2(3)*(school.').*repmat(ntf,1,R)+theta2(4)*(plincome.').*repmat(ntf,1,R);
% theta2(end-1)*ppi'.*repmat(price,1,R),theta2(2)*(dhnet.').*repmat(ntf,1,R),theta2(4)*(female.').*repmat(nss,1,R)

while dst > tol
        U = repmat(delta,1,R) + mu;
        eU = exp(U); %jm-by-R
        
        denom0 = 1 +sum(eU,1);
        denom = repmat(denom0,nj,1);
        
        P = mean(eU./denom,2);
        delta_new = delta + logS - log(P);
        dst = max(abs(delta - delta_new));
        delta = delta_new;
end
    

% while dst > tol
%     U = (repmat(delta,1,R) + mu)./(1-rho);
%     eU = exp(U); %j-by-R
%     
%     nestmatt = reshape(nestmat,[nj,1,nj]);
%     nestmati = repmat(nestmatt,1,R,1);
%     
%     nsteU = nestmati.*eU;
%     nsteUsum = sum(nsteU,3); % j-by-R
%     nsteUsum1 = nstidx.*nsteUsum;
%     
%     % Iig: j-by-R
%     Iig = (1-rho).*log(nsteUsum);
%     
%     % Ii: j-by-R
%     Iivec = log(1+sum(nsteUsum1,1)); % 1-by-R
%     Ii = repmat(Iivec,nj,1); % j-by-R
%    
%     num = eU.*exp(Iig);
%     denom = exp(Iig./(1-rho)).*exp(Ii);
%     
%     Pi = num./denom;
%     
%     P = mean(Pi,2);
%     logP = log(P);
%     logS = log(S);
%     delta_new = delta + logS - (1-rho)*log(P);
%     dst = max(abs(delta - delta_new));
%     delta = delta_new;
% end


end




