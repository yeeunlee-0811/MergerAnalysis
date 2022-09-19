function li = pyihat_wog(theta,inp)

nu = inp.nu;  % nu: individual Unovservables (n by R)
n = inp.n;
R = inp.R;

% individual vars
jdum = inp.jdum;
payer = inp.payer;
school = inp.school;
djob = inp.djob;
age = inp.age;
hhldsiz = inp.hhldsiz;
ichngnet3 = inp.ichngnet3;
female = inp.female;
dmar = inp.dmar;
income = inp.income;
linc = inp.linc;
dhnet = inp.dhnet;
dtviv = inp.dtviv;
dwviv = inp.dwviv;

% product vars
nj0 = inp.nj0;
pkmat_wog = inp.pkmat_wog; % nj0 x npk_wog (=3)
npk_wog = inp.npk_wog; % npk_wog = 3
netflix = inp.netflix;
wavve = inp.wavve;
tving = inp.tving;

thpk = theta(1:npk_wog,1); % npk_wog = 3
thnu = theta(npk_wog+1:npk_wog+3,1); % 6  
thnu1 = thnu(1,1);
thnu2 = thnu(2,1);
thnu3 = thnu(3,1);

this = theta(npk_wog+4:end,1); % 7~
th1 = this(1,1);
th2 = this(2,1);
th3 = this(3,1); 
th4 = this(4,1); 
th5 = this(5,1); 
th6 = this(6,1); 
th7 = this(7,1); 
th8 = this(8,1); 
th9 = this(9,1); 
th10 = this(10,1); 
th11 = this(11,1); 
th12 = this(12,1); 
th13 = this(13,1); 
th14 = this(14,1); 
th15 = this(15,1); 
% th16 = this(16,1); 
% th17 = this(17,1); 
% th18 = this(18,1); 
% th19 = this(19,1); 
% th20 = this(20,1); 
% th21 = this(21,1); 
% th22 = this(22,1); 
% th23 = this(23,1); 
% th24 = this(24,1); 
% th25 = this(25,1); 
% th26 = this(26,1); 
% th27 = this(27,1); 


pkmati0 = reshape(pkmat_wog,[nj0,1,npk_wog]);
pkmati = repmat(pkmati0,1,n); % nj0 x n x npk_wog

payermat = repmat(payer',1,1,R); % 1x n x R
dhnetmat = repmat(dhnet',1,1,R); 
dtvivmat = repmat(dtviv',1,1,R);
dwvivmat = repmat(dwviv',1,1,R);


agemat = repmat(age',1,1,R);
schoolmat = repmat(school',1,1,R);
djobmat = repmat(djob',1,1,R);
femalemat = repmat(female',1,1,R);
dmarmat = repmat(dmar',1,1,R);
incomemat = repmat(income',1,1,R);
lincmat = repmat(linc',1,1,R);
ihhmat = repmat(ichngnet3',1,1,R);
nu = reshape(nu,[1,n,R]);
hhldmat = repmat(hhldsiz',1,1,R);



mu = 0;
for kk = 1:npk_wog
    pki = pkmati(:,:,kk); % nj0 x n x npk_wog
    mu = mu + pki*thpk(kk,1); % nj0 x n
    munu0 = repmat(mu,1,1,R); % nj0 x n x R
end
    
    munu = munu0 ... 
        + thnu1*netflix.*nu +thnu2*wavve.*nu +thnu3*tving.*nu... 
        + th1*wavve.*payermat + th2*wavve.*dwvivmat + th3*wavve.*agemat + th4*wavve.*schoolmat + th5*wavve.*djobmat ... 
        + th6*netflix.*payermat + th7*netflix.*dhnetmat +th8*netflix.*agemat + th9*netflix.*schoolmat + th10*netflix.*djobmat...
        + th11*tving.*payermat + th12*tving.*dtvivmat + th13*tving.*agemat + th14*tving.*schoolmat + th15*tving.*djobmat;



expU = exp(munu); % nj0 x n x R
denom0 = sum(expU,1); % 1 x n x R
denom = repmat(denom0,nj0,1); % nj0 x n x R
phat0i = expU./denom; % nj0 x n x R
tjdum = jdum'; % nj0 x n
phati = phat0i.*tjdum; % nj0 x n x R

li0 = mean(phati,3); % nj0 x n x 1 

li = sum(li0,1); % 1 x n x 1


end