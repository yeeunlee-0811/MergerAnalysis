function gi = gi(theta,inp)

Z = inp.Z;
nj = inp.nj;

m = size(Z,2);

u = residualFun(theta,inp);

gi = zeros(nj,m);

for i =1:m
    gi(:,i)= Z(:,i).*u;
end

end
