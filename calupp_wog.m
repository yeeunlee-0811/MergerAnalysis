function upp = calupp_wog(theta,inp)

costvec = inp.costvec; % 3 x 1 , Order: net, wavve, tving (premerger cost)
S = inp.Sfora;
j = inp.pch;
k = inp.psa;
alpha = inp.alpha;
price = inp.price;
e1 = inp.e1;
e2 = inp.e2;

q0 = recalq(S); % 4 x 1
q = q0(2:end,1); % 3 x 1

alphaj = alpha(j,1);

pk = price(k,1);

cj = costvec(j,1);
ck = costvec(k,1);

dqdd0= diffdelta_wog(theta,inp); % 4 x 3
dqdd = dqdd0(2:end,:); % 3 x 3
dqkddj = dqdd(k,j); % p = p_j & q = q_k, scalar
dqjddj = dqdd(j,j); % p = p_j & q = q_j, scalar
upp0 = -(pk-(1-e1)*ck)*dqkddj;
upp1 = dqjddj;
upp2 = - e1*(cj);

upp = upp0/upp1+upp2;
end