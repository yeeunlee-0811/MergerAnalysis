function [fnval] = diffdelta(theta,inp)
% calculating dq/ddelta
% 1. parameter value changes (first three of thpk)
% 2. individual pick one which gives the best utility
% 3. How do we handle the nu? Pick one which gives the best average utility

preshare = sharefn(theta,inp); % nj0 x 1
s0 = recalq(preshare); % 4 x 1

fnval0 = zeros(4,3); 
s1mat = zeros(4,3);
thetanew = theta;

for k = 1:3
   thetanew(k,1) = theta(k,1) + 0.05;
   postshare = sharefn(thetanew,inp);
   s1 = recalq(postshare); % 4 x 1
   s1mat(:,k) = s1;
   fnval0(:,k) = (s1 - s0)./0.05;
   thetanew = theta;
end

fnval = fnval0;
end


