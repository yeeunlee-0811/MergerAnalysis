function [fnval] = calalpha(alpha,theta0,price,inp)

S = inp.Sfora; % real share or predicted share (nj0 x 1)
costvec = inp.costvec; % 3 x 1

s4 = recalq(S); % 4 x 1
s = s4(2:end,1); % 3 x 1

dqddmat0 = diffdelta(theta0,inp); % 4 x 3
dqddmat = dqddmat0(2:end,:);
own0 = inp.own0;

fnval = s - alpha.*diag(own0.*dqddmat).*(price-costvec);
end

