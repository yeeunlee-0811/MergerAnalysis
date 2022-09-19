clear 

%% Data Import
data0 =  readtable('data.csv');
data=table2array(data0);

%%
pcharac0 = readtable('pcharac.xlsx');
xpk = table2array(pcharac0);


%% Restricting Market

% Market Index = 1 market, yr = 2021
yr = data0.yr;
iyr = yr ==21;
data21 = data(iyr,:);

%% Generating variables 
n = size(data21,1);
kt = size(data21,2);
varnames = data0.Properties.VariableNames;

%%
Xpk = xpk(2:end,:);

ncpt = size(xpk,2);
nj0=size(xpk,1);
nj = nj0-1;

pkvarnames = pcharac0.Properties.VariableNames;

%% Assigning variable names

for i = 1:kt
    eval([varnames{i} '=  data21(:,i)']);
end

for i = 1:ncpt
    eval([pkvarnames{i} '=  Xpk(:,i)']);
end

%% Variables

% Market share of goods
inp.nj0 = nj0;
inp.nj = nj;

% goods dummy
% dgk: indicates goods that without netflix but consuming kotts
dgk = data21(:,[140:143, 167, 145:150, 163:166]);

% dgnk: indicates goods that with netflix and kotts
idn = ottnet == 1;
dn = double(idn); % netflix user 
dn0 = ones(n,1)-dn;
dgn1k1 = dn.*dgk; % netflix + kotts
dgn0k1 = dn0.*dgk;

% dgn: ppl who are using only Netflix
in1k1 = sum(dgn1k1,2)> 0;
in1k1 = double(in1k1);
dgn = dn-in1k1;

% dy1 : people who are using only ytub
iy1 = nott ==1 & ott1 == 10;
dy1 = double(iy1);

idytub1 = dytub == 1;
dytub1 = double(idytub1);
dytub0 = ones(n,1)-dytub1;

%y0n0k0
iy0n0k0 = sum([dy1,dgk,dgn],2)==0;
dy0n0k0 = double(iy0n0k0);

mynk = [dy0n0k0,dytub0.*dgn0k1,dytub0.*dgn,dytub0.*dgn1k1,dy1,dytub1.*dgn0k1,dytub1.*dgn,dytub1.*dgn1k1];


nj = nj0-1; %outside option remove
inp.nj = nj;

Slevel = sum(mynk,1);
S0 = Slevel./n;
S1 = S0';
S2 = S1(2:end,1); 

%% Removing zero share products
iS1 = S2>0;
S = S2(iS1,1);
inp.S = S;

nj = size(S,1);
inp.nj = nj;

%% Nest
% Total 8 nests
nestvec0 = [ones(5,1);2*ones(10,1);3*ones(6,1);4*ones(10,1);5*ones(6,1);6*ones(10,1);7*ones(6,1);8*ones(10,1)];
nestvec = nestvec0(iS1,1);

inp.nestvec = nestvec;

nestmat = zeros(nj,nj);

for i = 1:nj
    nst = nestvec(i,1);
    inst = nestvec == nst;
    dinst = double(inst);
    nestmat(i,:) = dinst.';
end

inp.nestmat = nestmat;% nj-by-nj

nstidx = zeros(nj,1);
nstidx([1,6,13,18,25,31,41],1) = 1;

inp.nstidx = nstidx;

%% RC unobservable


R = n; % number of random draws
inp.R = R;
% Random coefficient only applied to price
%Krc = 1;
%inp.Krc = Krc;


% different draws for each year, and each product characteristic.
rng(10) 
nu0 = randn(n,3); %Market level draws
%nu = ones(nj,1).*nu0;

inp.nu =nu0;

%% Product characteristics

xxpk=Xpk(iS1,:);
for i = 1:ncpt
    eval([pkvarnames{i} '=  xxpk(:,i)']);
end


dkss = ones(nj,1);

ik120 = k1 == 0 & k2 == 0;
dkss(ik120,:) =0;

inp.dkss = dkss;

dmm = dm2+dm3;
ddm = dd2+dd3;
ntf = netflix;
dkk = koreanss;
dytb1 = dytb -dynk-dyn-dyk;
dntf1 = ntf - dynk-dyn-dnk;
dkss1 = dkss -dynk-dyk-dnk;

dm = dm1 + dm2 + dm3;
dm23 = dm2+dm3;
db = db1+ db2;

inp.dm = dm;
inp.db = db;

jdum = [dytb,ntf,dkss];
prodcharac0 = [ones(nj,1),price,dtv,dwv,dktv,dmb,dns,dytb,dynk,dnk,dkk];
%dmb,dns
inp.prodcharac0 = prodcharac0;

inp.price = price;
inp.avod = avod;
inp.tvod = tvod;
inp.nss = nss;
inp.fdrama = fdrama;
inp.broadcast = broadcast;
inp.dkss = dkss;
inp.ntf = ntf;
inp.dm1 = dm1;
inp.dm2 = dm2;
inp.dm3 = dm3;
inp.J = J;
inp.davod = davod;
inp.dsvod = dsvod;
inp.dtvod = dtvod;
inp.dfdrama = dfdrama;
inp.db1 = db1;
inp.db2 = db2;

inp.dd1 = dd1;
inp.dd2 = dd2;
inp.dd3 = dd3;

inp.dytb = dytb;

%% check valid prodcharac
Zx = prodcharac0;
zp = Zx(:,1:2);
for i=3:size(Zx,2)
    zi = Zx(:,i);
    reg = fitlm(zp,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        i
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi);
        end
        zp = [zp zi];
        % disp(i)
    end
end

%%


prodcharac = prodcharac0(:,2:end);
inp.prodcharac=prodcharac;

%% Z

Zx = [prodcharac0];
zp = Zx(:,1:2);
for i=3:size(Zx,2)
    zi = Zx(:,i);
    reg = fitlm(zp,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        i
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi);
        end
        zp = [zp zi];
        % disp(i)
    end
end
%%
%%
maxnss = max(nss); % = 4

zp0 = [prodcharac];

% average price of a group of the same notts
instr0 = zp0; 
ninst = size(instr0,2);

Z0 = zeros(nj,ninst,2); % own + in-nest
for i= 1: nj
    ssn = nss(i,1);
    sames = nss == ssn;
    nsames = nss ~= ssn;

    for j =1:ninst
    tinst = instr0(:,j);
    Z0(i,j,1) =tinst(i);
    Z0(i,j,2) = sum(tinst(nsames));
    end
end

Z10 = reshape(Z0,nj,ninst*2);
Z1 = Z10(:,2:end);
%%
Z_0 = [ones(nj,1),Z1] ;

z0_D = Z_0(:,1:2);
for i=3:size(Z_0,2)
    zi = Z_0(:,i);
    reg = fitlm(z0_D,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        i
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi);
        end
        z0_D = [z0_D zi];
        % disp(i)
    end
end

Zfin = z0_D;
%inp.Z_blp=Z_blp;

%%
Zxx = [Zfin];
zp = Zxx(:,1:2);
for i=3:size(Zxx,2)
    zi = Zxx(:,i);
    reg = fitlm(zp,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        i
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi);
        end
        zp = [zp zi];
        % disp(i)
    end
end
%%

Z_1 = [Zfin,Zfin(:,3).*Zfin(:,2),Zfin(:,3).*Zfin(:,4)] ;
%Zfin(:,5).*Zfin(:,6),Zfin(:,6).*Zfin(:,7),Zfin(:,7).*Zfin(:,8),Zfin(:,3).*Zfin(:,4)

z1_D = Z_1(:,1:2);
for i=3:size(Z_1,2)
    zi = Z_1(:,i);
    reg = fitlm(z1_D,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        i
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi);
        end
        z1_D = [z1_D zi];
        % disp(i)
    end
end

Zfin1 = z1_D;
%inp.Z_blp=Z_blp;

Z = Zfin1;
inp.Z= Z;

%%

payer = ones(n,1);
ip0 = pmtss == 0;
payer(ip0,1) = 0;

ipin = pi<0;
ipip = pi>0;
dipin = zeros(n,1)-ones(n,1).*double(ipin);
dipip = zeros(n,1)+ones(n,1).*double(ipip);
dpi = dipin+dipip;

inp.payer = payer;
%%
inp.school = school;
inp.female = female;
inp.dmar = dmar;
inp.djob = djob;
inp.dytub = dytub;
inp.ij = ij;
inp.n = n;
inp.payer = payer;
inp.pi = pi;

inet = zeros(n,1);
iinet = 1 == ottnet;
diinet = double(iinet);
netflix = inet + diinet;
inp.ihhnet = ichngnet3;

ikss = zeros(n,1);
iikss = 1 == dkott;
diikss = double(iikss);
koreanss = ikss + diikss;

inp.netflix = netflix;
inp.koreanss = koreanss;
inp.youtube = dytub1;

indivd = [dpi,payer.*income,hhldsiz,age,dhnet];
inp.indivd = indivd;
Kindiv = size(indivd,2)+1;
inp.Kindiv = Kindiv;

inanmat = zeros(n, Kindiv-1);

for k = 1:Kindiv-1
   inan = isnan(indivd(:,k));
   dinan = double(inan);
    
   inanmat(:,k) = dinan;
end

sum(inanmat,1)

%%

W = (Z'*Z)/nj; % initial weighting matrix
inp.W = W;
inp.conc = true; % "concentrate out" linear parameters


residual_fn = @(theta) residualFun(theta,inp); 
obj = @(beta) GMM_objective(beta,Z,W,residual_fn);

theta0 = [0.0920896170988548;0.509704901478813;0.5];
%rand(2,1); % starting values; randomly drawn from (0,1)
% [0.145758053583976;0.762083183019301] : starting value of pincome and
% agenss

% %% Run fminuncf
% options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
%     'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
% theta2_0 = fminunc(obj,theta0,options); 

options = optimset('Display','iter','MaxFunEvals',20000);
theta2_0 = fminsearch(obj,theta0,options);
[~,theta1_0]=residualFun(theta2_0,inp); % theta1 = parameters in the mean util, # 1+#prodcharac+premium
theta_0 = [theta1_0;theta2_0]; 

%%
% options = optimset('Display','iter','MaxFunEvals',20000);
% theta2_0 = fminsearch(obj,theta0,options);

%%
inp.conc =0;
[var] = sefn0(theta_0,inp);
se_ = sqrt(diag(var));
tstat = theta_0./se_;

%%
save('se.mat','se')
save('tstat.mat','tstat')
save('theta_0.mat','theta_0')

%% What if Netflix price increases 1$
% Calculating P with the parameter estimate
price1 = price;
intf = ntf == 1;
price(intf,:) = price(intf,:)+1;


P = calcs(price1, theta_0,inp);

%% 
spj = [dytb,ntf,k1,k2];
