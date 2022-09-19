function s = recalq(q)
% q: nj0 x 1 vector
% s: 4 x 1

q0 = q(1,1);
q1 = q(2,1);
q2 = q(3,1);
q3 = q(4,1);
q4 = q(5,1);
q5 = q(6,1);
q6 = q(7,1);
q7 = q(8,1);

s = zeros(4,1);
s(1,1) = q0;
s(2,1) = q1 + q5 + q6 + q7;
s(3,1) = q2 + q4 + q6 + q7;
s(4,1) = q3 + q4 + q5 + q7;

end


