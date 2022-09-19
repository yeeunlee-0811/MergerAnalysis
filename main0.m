clear 

%% Data Import
data0 =  readtable('dtrev.csv');
data = table2array(data0);

%% Restricting year: 2020 and 2021
yr = data0.yr;
iyr = yr == 21;
datar = data(iyr,:);

%% Generating individual variables 
n = size(datar,1);
nikt = size(datar,2);
varnames = data0.Properties.VariableNames;

% Assigning individual variable names

for i = 1:nikt
    eval([varnames{i} '=  datar(:,i)']);
end

%%
dothers = zeros(n,1);
didothers = nott >0;
ddidothers = double(didothers);
dothers(didothers,1) = 1 ;

servp = zeros(n,1);
didservp = serviceprovider == 1;
servp(didservp,1) = 1;

conti = zeros(n,1);
didconti = ottcontents121 == 1 | ottcontents121 == 3;
conti(didconti,1) = 1;

direl = zeros(n,1);
didirel = irel == 1;
direl(didirel,1) = 1;

%%
inp.conti = conti;
inp.dothers = dothers;
inp.n = n;
inp.payer = payer;
inp.school = school;
inp.female = female;
inp.dmar = dmar;
inp.djob = djob;
inp.income = income;
inp.age = age;
inp.hhldsiz = hhldsiz;
inp.nott = nott;
inp.ichngnet3 = ichngnet3;
inp.dhnet = dhnet;
inp.dtviv = dtviv;
inp.dwviv = dwviv;
inp.linc = log(income);
inp.servp = servp;
inp.direl = direl;

%% Generating product variables

pcharac0 = readtable('pcharac_new.xlsx');
xpk = table2array(pcharac0);

pkmat0 = xpk;
pkmat = pkmat0(:,2:end);

npk = size(pkmat,2); % number of products - combinations
nj0 = size(pkmat0,1);
nj = nj0-1;


pkvarnames = pcharac0.Properties.VariableNames;

for i = 1:npk
    j = i+1;
    eval([pkvarnames{j} '=  pkmat(:,i)']);
end

%%
inp.nj0 = nj0;
inp.nj = nj;
inp.pkmat=pkmat;
inp.npk = npk;
inp.netflix = netflix;
inp.wavve = wavve;
inp.tving = tving;
inp.wvtv = wvtv;
inp.nettv = nettv;
inp.netwv = netwv;
inp.nwt = nwt;

%% Market Share 

jdum = zeros(n,nj0);
for ij = 1:nj0
    ij0 = ij - 1;
    iji = catvar == ij0;
    diji = double(iji);
    jdum(:,ij) = diji;
end

Slevel = sum(jdum,1);
S0 = Slevel./n;
S = S0';

%% 
inp.S = S;
inp.jdum = jdum;

%% Simulated Draws

R = 50; % R: The number of draws
inp.R = R;

% Individual Level random draws
rng(10) 
nu = randn(n,R); % n x R
save('nu.mat','nu')

%%
load('nu.mat');
inp.nu = nu;


%% MLE

pcalalpha = [9.26; 8.41; 8.41];
inp.price = pcalalpha;

llobj = @(theta)loglikelihood(theta,inp);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',1e10000);

% N. of parameters? 
np = 34;
inp.np = np;

param0 = 0.5*ones(np,1);
[mleparam,logval] = fminsearch(llobj,param0,options);

%% Getting variance
%mleparam = [-7.099; -4.844;-2.479;3.392;-2.927;3.457;-2.882;2.062;-2.059;1.393;2.939;6.511;-0.135;-0.892;-1.042;4.252;4.358;-0.832;1.133;0.351;2.085;4.036;-0.732;0.128;0.981];

vars_matrix = mlevars(mleparam,inp);
se = sqrt(diag(vars_matrix));



%% Test sharefn

costvec = [7.35; 6.12; 6.12]; % net, wav, tv
costavg = mean(costvec);
costvec_avg = costavg*ones(3,1);

% adjust costvec input
inp.costvec = costvec;

% test sharefn and compare predicted share (preS) and actual share(S)
preS = sharefn(mleparam,inp);
pres4 = recalq(preS);
s4 = recalq(S);
%% calculate alpha
own0 = eye(3,3);
inp.own0 = eye(3,3);
inp.Sfora = preS; % real share or predicted share (nj0 x 1)
calalphafn = @(alpha) calalpha(alpha,mleparam,pcalalpha,inp);

alpha0 = 0.5*ones(3,1); 
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-12,...
    'MaxIterations',25);
[gotalpha,fval0] = fsolve(calalphafn,alpha0,options); % starting value = observed prices

%% .........................................using the parameter alpha, do the merger simulation - 1...............................................................
% finding price that makes the three FOCs are equal to zero
e1 = 0.3;
inp.e1 = e1;
e2 = 0.3;
inp.e2 = e2;

costvecpm = costvec.*[1;1-e1;1-e2];
inp.costvecpm = costvecpm; %%%%% check what costvector you using

own = [1, 0, 0; 0, 1, 1; 0, 1, 1];
inp.own = own;
alpha = gotalpha; % 3 x 1
inp.alpha = alpha;
deltabar = mleparam(1:3,1) + alpha.*pcalalpha;
inp.deltabar = deltabar;

%%
msfn = @(price) profitfoc(price, mleparam, inp);


price0 =[4.14391101185198;9.50702045750427;2.88767283368712];
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-100,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-3,...
    'MaxIterations',2500000000000000);
[mergerp,fval1] = fsolve(msfn,price0,options); % starting value = observed prices

%%
[~,pmshare] = profitfoc(mergerp, mleparam, inp);
pms4 = recalq(pmshare);
pms43 = pms4(2:end,1);
profitpm1 = (mergerp-costvecpm).*pms43;
pres43 = pres4(2:end,1);
profitprem = (pcalalpha-costvec).*pres43;

%% using the parameter alpha, calculate the UPP
% net, wav, tv
pch = 3;
inp.pch = pch;
psa = 2;
inp.psa = psa;

upp32 = calupp(mleparam,inp);
%%
pch = 2;
inp.pch = pch;
psa = 3;
inp.psa = psa;


upp23 = calupp(mleparam,inp);

%% NO-Bundle

e1 = 0.3;
inp.e1 = e1;
e2 = 0.3;
inp.e2 = e2;
costvecpm = costvec.*[1;1-e1;1-e2];
inp.costvecpm = costvecpm; %%%%% check what costvector you using

pcharac0_nob = readtable('pcharac_nob.xlsx');
xpk_nob = table2array(pcharac0_nob);

pkmat0_nob = xpk_nob;
pkmat_nob = pkmat0_nob(:,2:end);

npk_nob = size(pkmat_nob,2); % number of products - combinations
nj0_nob = size(pkmat0_nob,1);
nj_nob = nj0_nob-1;


pkvarnames_nob = pcharac0_nob.Properties.VariableNames;

for i = 1:npk_nob
    j = i+1;
    eval([pkvarnames_nob{j} '=  pkmat_nob(:,i)']);
end

%%
inp.nj0 = nj0_nob;
inp.nj = nj_nob;
inp.pkmat=pkmat_nob;
inp.npk = npk_nob;
inp.netflix = netflix_nob;
inp.wavve = wavve_nob;
inp.tving = tving_nob;
inp.wvtv = wvtv_nob;
inp.nettv = nettv_nob;
inp.netwv = netwv_nob;
inp.nwt = nwt_nob;

%%
% % Adjust variables
inp.own_nob = eye(3,3);

inp.costvec = costvecpm;

% % 
msfn_nob = @(price) profitfoc_nob(price, mleparam, inp);

price0 = [11.3121878399901;2.91965568472165;9.57493975286354];

inp.price = pcalalpha;
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-100,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-4,...
    'MaxIterations',2500000000000000);
[mergerp_nob,fval_nob] = fsolve(msfn_nob,price0,options); % starting value = observed prices

[~,pmshare_nob] = profitfoc_nob(mergerp_nob, mleparam, inp);

pms4_nob = recalq_nob(pmshare_nob);
pms43_nob = pms4_nob(2:end,1);
profitpm1_nob = (mergerp_nob-costvecpm).*pms43_nob;

%% Simulation when merged firm sells only one good

%% Combined one
e1 = 0;
inp.e1 = e1;
e2 = 0;
inp.e2 = e2;
costvecpm = costvec.*[1;1-e1;1-e2];
inp.costvecpm = costvecpm; %%%%% check what costvector you using

pcharac0_cg = readtable('pcharac_cg.xlsx');
xpk_cg = table2array(pcharac0_cg);

pkmat0_cg = xpk_cg;
pkmat_cg = pkmat0_cg(:,2:end);

npk_cg = size(pkmat_cg,2); % number of products - combinations
nj0_cg = size(pkmat0_cg,1);
nj_cg = nj0_cg-1;


pkvarnames_cg = pcharac0_cg.Properties.VariableNames;

for i = 1:npk_cg
    j = i+1;
    eval([pkvarnames_cg{j} '=  pkmat_cg(:,i)']);
end

%%
inp.nj0 = nj0_cg;
inp.nj = nj_cg;
inp.pkmat=pkmat_cg;
inp.npk = npk_cg;
inp.netflix = netflix_cg;
inp.wavve = wavve_cg;
inp.tving = tving_cg;
inp.wvtv = wvtv_cg;
inp.nettv = nettv_cg;
inp.netwv = netwv_cg;
inp.nwt = nwt_cg;

%%
% % Adjust variables
inp.own_cg = eye(2,2);

cost1 = costvecpm(1,1);
cost2 = costvecpm(2,1);
cost3 = costvecpm(3,1);
cost23 = (cost2+cost3);

costvec_cg = [cost1;cost23];
inp.costvec = costvec_cg;

% % 
ms1fn = @(price) profitfocs_two(price, mleparam, inp);

p1 = pcalalpha(1,1);
p2 = mean(pcalalpha(2:end,1));

price0 = abs(randn(2,1))*10;
% [7.86570098948843;8.28464283069710] - delta2 change
% [3.61202571292274;3.13809198805179] - delta3 change

inp.price = pcalalpha;
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-100,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-4,...
    'MaxIterations',2500000000000000);
[mergerp1,fval1] = fsolve(ms1fn,price0,options); % starting value = observed prices

[~,pmshare1] = profitfocs_two(mergerp1, mleparam, inp);

pms4_1 = recalq_two(pmshare1);
pms43_1 = pms4_1(2:end,1);
profitpm1_1 = (mergerp1-costvec_cg).*pms43_1;

%% One good - tving
pcharac0_pmt = readtable('pcharac_pmt.xlsx');
xpk_pmt = table2array(pcharac0_pmt);

pkmat0_pmt = xpk_pmt;
pkmat_pmt = pkmat0_pmt(:,2:end);

npk_pmt = size(pkmat_pmt,2); % number of products - combinations
nj0_pmt = size(pkmat0_pmt,1);
nj_pmt = nj0_pmt-1;


pkvarnames_pmt = pcharac0_pmt.Properties.VariableNames;

for i = 1:npk_pmt
    j = i+1;
    eval([pkvarnames_pmt{j} '=  pkmat_pmt(:,i)']);
end

%%
inp.nj0 = nj0_pmt;
inp.nj = nj_pmt;
inp.pkmat=pkmat_pmt;
inp.npk = npk_pmt;
inp.netflix = netflix_pmt;
inp.wavve = wavve_pmt;
inp.tving = tving_pmt;
inp.wvtv = wvtv_pmt;
inp.nettv = nettv_pmt;
inp.netwv = netwv_pmt;
inp.nwt = nwt_pmt;

%% Tving
% % Adjust variables
inp.own_pmt = eye(2,2);

e1 = 0;
inp.e1 = e1;
e2 = 0;
inp.e2 = e2;
costvecpm = costvec.*[1;1-e1;1-e2];
inp.costvecpm = costvecpm; %%%%% check what costvector you using

cost1 = costvecpm(1,1);
cost2 = costvecpm(2,1);
cost3 = costvecpm(3,1);
cost23 = (cost2+cost3)/2;

costvec_pmt = [cost1;cost3];
inp.costvec = costvec_pmt;

p1 = pcalalpha(1,1);
p2 = mean(pcalalpha(2:end,1));

% % 
ms2pmtfn = @(price) profitfocs_two(price, mleparam, inp);

price0 = abs(randn(2,1))*10;

inp.price = pcalalpha;
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-100,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-4,...
    'MaxIterations',2500000000000000);
[mergerp2pmt,fval2pmt] = fsolve(ms2pmtfn,price0,options); % starting value = observed prices

[~,pmshare2pmt] = profitfocs_two(mergerp2pmt, mleparam, inp);

pms4_2pmt = recalq_two(pmshare2pmt);
pms43_2pmt = pms4_2pmt(2:end,1);
profitpm2_pmt = (mergerp2pmt-costvec_pmt).*pms43_2pmt;


%% One good - wavve
pcharac0_pmw = readtable('pcharac_pmw.xlsx');
xpk_pmw = table2array(pcharac0_pmw);

pkmat0_pmw = xpk_pmw;
pkmat_pmw = pkmat0_pmw(:,2:end);

npk_pmw = size(pkmat_pmw,2); % number of products - combinations
nj0_pmw = size(pkmat0_pmw,1);
nj_pmw = nj0_pmw-1;


pkvarnames_pmw = pcharac0_pmw.Properties.VariableNames;

for i = 1:npk_pmw
    j = i+1;
    eval([pkvarnames_pmw{j} '=  pkmat_pmw(:,i)']);
end

%%
inp.nj0 = nj0_pmw;
inp.nj = nj_pmw;
inp.pkmat=pkmat_pmw;
inp.npk = npk_pmw;
inp.netflix = netflix_wv;
inp.wavve = wavve_wv;
inp.tving = tving_wv;
inp.wvtv = wvtv_wv;
inp.nettv = nettv_wv;
inp.netwv = netwv_wv;
inp.nwt = nwt_wv;

%% Wavve
% % Adjust variables
inp.own_pwv = eye(2,2);
e1 = 0.3;
inp.e1 = e1;
e2 = 0.3;
inp.e2 = e2;
costvecpm = costvec.*[1;1-e1;1-e2];
inp.costvecpm = costvecpm; %%%%% check

cost1 = costvecpm(1,1);
cost2 = costvecpm(2,1);
cost3 = costvecpm(3,1);
cost23 = (cost2+cost3)/2;

costvec_pmw = [cost1;cost2];
inp.costvec = costvec_pmw;

% % 
ms2pmwfn = @(price) profitfocs_two(price, mleparam, inp);

price0 = abs(randn(2,1))*10;

inp.price = pcalalpha;
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-100,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-4,...
    'MaxIterations',2500000000000000);
[mergerp2pmw,fval2pmw] = fsolve(ms2pmwfn,price0,options); % starting value = observed prices

[~,pmshare2pmw] = profitfocs_two(mergerp2pmw, mleparam, inp);

pms4_2pmw = recalq_two(pmshare2pmw);
pms43_2pmw = pms4_2pmw(2:end,1);
profitpm2_pmw = (mergerp2pmw-costvec_pmw).*pms43_2pmw;


%% .......................Without Gamma terms...................................................................................
%% product characteristics _ wog

pkmat_wog = pkmat(:,1:3);
npk_wog = size(pkmat_wog,2); % number of products - combinations

inp.pkmat_wog = pkmat_wog;
inp.npk_wog = npk_wog;

% % Use only Netlfix, Wavve and Tving 
% inp.netflix = netflix;
% inp.wavve = wavve;
% inp.tving = tving;

% % Do not use these bundle good terms
% inp.wvtv = wvtv;
% inp.nettv = nettv;
% inp.netwv = netwv;
% inp.nwt = nwt;

%% MLE _ WO gamma

llobj_wog = @(theta)loglikelihood_wog(theta,inp);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',1e10000);

% N. of parameters? 
npwg = 21;
inp.npwg = npwg;

param0_wog = 0.5*ones(npwg,1);
[mleparam_wog,logval] = fminsearch(llobj_wog,param0_wog,options);
%%
mleparam_wog = [-13.9002427653501;-3.21315699553498;-3.12357073798712;-2.30051630156088;9.84383888751971;0.270571348106650;5.78415779101309;11.5047602544210;-3.77966692829957;-0.667819840170335;4.53063687058438;1.47074733428552;5.91385061499166;-0.439022712832211;2.29409651589224;-1.29399089225768;0.663992275007238;3.44793804457806;-0.432006535317514;0.135386122066149;0.602415732664724];

% Getting Variance
vars_matrix_wog = mlevars_wog(mleparam_wog,inp);
se_wog = sqrt(diag(vars_matrix_wog));

%% calculate alpha

costvec = [7.35; 6.12; 6.12]; % net, wav, tv
costvec_avg = 6.27*ones(3,1);

% adjust costvec input
inp.costvec = costvec;

% test sharefn and compare predicted share (preS) and actual share(S)
preS_wog = sharefn_wog(mleparam_wog,inp);
pres4_wog = recalq(preS_wog);
s4 = recalq(S);
%%
pcalalpha = [9.26; 8.41; 8.41];
inp.price = pcalalpha;
own0 = eye(3,3);
inp.own0 = eye(3,3);
inp.Sfora = preS_wog; % real share or predicted share (nj0 x 1)
calalphafn_wog = @(alpha) calalpha_wog(alpha,mleparam_wog,pcalalpha,inp);

alpha0 = ones(3,1); 
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-12,...
    'MaxIterations',25);
[gotalpha_wog,fval0] = fsolve(calalphafn_wog,alpha0,options); % starting value = observed prices

% gotalpha_wog = [1.08221024258760;0.503968253968252;0.371443624868283];
%% using the parameter alpha, do the merger simulation
% finding price that makes the three FOCs are equal to zero
e1 = 0.3;
inp.e1 = e1;
e2 = 0.3;
inp.e2 = e2;

costvecpm = costvec.*[1;1-e1;1-e2];
inp.costvecpm = costvecpm; %%%%% check what costvector you using

own = [1, 0, 0; 0, 1, 1; 0, 1, 1];
inp.own = own;
alpha_wog = gotalpha_wog; % 3 x 1
inp.alpha_wog = alpha_wog;
deltabar_wog = mleparam_wog(1:3,1) + alpha_wog.*pcalalpha;
inp.deltabar_wog = deltabar_wog;

%%

msfn_wog = @(price) profitfoc_wog(price, mleparam_wog, inp);


price0 = pcalalpha; 
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-100,...
    'StepTolerance',1e-100,'FunctionTolerance',1e-4,...
    'MaxIterations',2500000000000000);
[mergerp_wog,fval1_wog] = fsolve(msfn_wog,price0,options); % starting value = observed prices


%% using the parameter alpha, calculate the UPP
% net, wav, tv
pch = 3;
inp.pch = pch;
psa = 2;
inp.psa = psa;

upp_wog32 = calupp_wog(mleparam_wog,inp);


%%
pch = 2;
inp.pch = pch;
psa = 3;
inp.psa = psa;

upp_wog23 = calupp_wog(mleparam_wog,inp);