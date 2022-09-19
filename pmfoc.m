function fnval = pmfoc(price, theta0, inp)

costvec = inp.costvecpm; % 3 x 1
deltabar = inp.deltabar;
alpha0 = inp.alpha;
dndp = inp.dndp; % 3 x 1 (price increase by 1) (share --> percentage)
dndq = inp.dndq; % 3 x 1 (dq: one percentage increase)

%% Reduce to two singleton goods - average out good 2 and good 3
alpha01 = alpha0(1,1);
alpha02 = alpha0(2,1);
alpha03 = alpha0(3,1);
alpha = [alpha01; (alpha02 + alpha03)./2];
av = - alpha; % 2 x 1


theta01 = theta0;
theta011 = theta01(1:3,1);
theta012 = theta01(4:end,1);
delta1 = theta011(1,1);
delta2 = theta011(2,1);
delta3 = theta011(3,1);
deltavec = [delta1; delta2+delta3];

deltabar = deltavec + alpha.*pcalalpha; % 2 x 1

theta11 = av.*price + deltabar; % now delta changes as price change
theta12 = theta012;
theta = [theta11;theta12]; %25 x 1

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
sec0 = av.*pricematt.*(own.*dqddt); % 3 x 3 
sec = sum(sec0,2); % 3 x 1 

% third term
thr = costvec.*dndp; % 3 x 1

% fourth term
costmat = costvec.*ones(3,3); % 3 x 3
costmatt = costmat';
dndqmat = dndq.*ones(3,3);
dndqmatt = dndqmat';
fou0 = av.*dndqmatt.*costmatt.*(own.*dqddt); % 3 x 3;
fou = sum(fou0,2); % 3 x 1

fnvalue = fir + sec - thr - fou;
% price
end
