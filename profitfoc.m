function [fnvalue,q] = profitfoc(price, theta0, inp)
%fnvalue set to be zero to get optimal price

costvec = inp.costvecpm; % 3 x 1
own = inp.own; % 3 x 3
deltabar = inp.deltabar;
alpha = inp.alpha;
n = 100;

av = - alpha;

theta = theta0;
theta(1:3,1) = av.*price + deltabar; % now delta changes as price change

% first term: share
q = sharefn(theta,inp); % nj0 x 1
sall = recalq(q); % 4 x 1
fir = sall(2:end,1); % 3 x 1

% dndp = (fir*n)./costvec; % 3 x 1 (price increase by 1) (share --> percentage)
% dndq = (price*n)./costvec; % 3 x 1 (dq: one percentage increase)

% second term
dqdd0 = diffdelta(theta,inp); % 4 x 3
dqdd1 = dqdd0(2:end,:); % 3 x 3
dqddt = dqdd1';
pricemat = price.*ones(3,3); % 3 x 3
pricematt = pricemat'; % 3 x 3
costmat = costvec.*ones(3,3); % 3 x 3
costmatt = costmat'; % 3 x 3
pminc = pricematt - costmatt;
sec0 = av.*pminc.*(own.*dqddt); % 3 x 3 
sec = sum(sec0,2); % 3 x 1 

fnvalue = fir + sec;
end











