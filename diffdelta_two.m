function dfnval = diffdelta_two(theta0,inp)
% diffdelta_two evaluaes dq/dd 
% the value needs to change when delta changes - initial point

theta_new = theta0;

q0 = sharefn_two(theta0,inp);
q03 = recalq_two(q0); % 3 x 1

dfnval = zeros(3,2);

k1 = inp.k1;
k2 = inp.k2;

theta_new(k1,1) = theta0(k1,1) + 0.5;
q1k1 = sharefn_two(theta_new,inp);
q13k1 = recalq_two(q1k1); % 3 x 1
dfnval(:,1) = (q13k1-q03)./0.5;

theta_new = theta0;

theta_new(k2,1) = theta0(k2,1) + 0.5;
q1k2 = sharefn_two(theta_new,inp);
q13k2 = recalq_two(q1k2); % 3 x 1
dfnval(:,2) = (q13k2-q03)./0.5;


end
