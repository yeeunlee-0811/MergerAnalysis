function [fnval,q0] = profitfocs_two(price,theta0,inp)

%% Adjust deltas to fit with the new setting

own = inp.own_cg;

delta2 = theta0(2,1);
delta3 = theta0(3,1);
delta23 = delta2 + delta3;

alpha = inp.alpha;
alpha1 = alpha(1,1);
alpha2 = alpha(2,1);
alpha3 = alpha(3,1);

av = -[alpha1;alpha3]; % change

price0 = inp.price;
price1 = price0(1,1);
price2 = price0(2,1);
price3 = price0(3,1);
pavg23 = mean(price0(2:end,1));

deltabar = inp.deltabar;
deltabar1 = deltabar(1,1);
deltabar2 = deltabar(2,1);
deltabar3 = deltabar(3,1);

p1 = price(1,1); % price is 2 x 1 vector
p2 = price(2,1);

dnew10 = -alpha1*p1 + deltabar1;
dnew20 = delta2;
% delta2;
% -alpha2*p2 + deltabar2;
dnew30 = -alpha3*p2 + deltabar3;
% delta3;
%  -alpha3*p2 + deltabar3;


%% Choose one of the merging firms
inp.k1 = 1;
inp.k2 = 3;

theta01 = theta0(1:3,1);
theta02 = theta0(4:end,1);

theta = [dnew10;dnew20;dnew30;theta02]; % the same number of thetas 

costvec = inp.costvec;
%% t1

q0 = sharefn_two(theta,inp); % sharefn_one generates total 4 shares: og + two singletons + one bundle
q1 = recalq_two(q0); % 3 x 1
q = q1(2:end,1); % 2 x 1

t1 = q; % 2 x 1

%% t2

dqdd0 = diffdelta_two(theta0,inp); % 3 x 2 
dqdd1 = dqdd0(2:end,:); % 2 x 2
sec0 = own.*dqdd1; % 2 x 2 
sec1 = diag(sec0); % 2 x 2 
sec2 = (price-costvec).*av.*sec1; % 2 x 1

t2 = sec2;


fnval = t1 + t2;

end





