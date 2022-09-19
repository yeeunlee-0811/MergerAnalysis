% this function calculates the 'xi'
function [res,theta1] = residualFun(theta,inp)

prodcharac = inp.prodcharac;
Z = inp.Z;
W = inp.W;
nj = inp.nj;

X = [ones(nj,1), prodcharac];
K = size(X,2);

if inp.conc % concentrating out linear parameters
    theta2 = theta;
else
    theta1 = theta(1:K); % linear parameters are put first.
    theta2 = theta(K+1:end);
end

delta = getdelta(theta2,inp); % only nonlinear parameters theta2 = (sigma) are inputs here.

Y = delta;

if inp.conc
    xzWz = X'*Z*(W\Z');
    A = xzWz*X;
    B = xzWz*Y;
    invA = inv(A);
    theta1 = A\B;
end
res = Y - X*theta1;

end