function q1 = recalq_two(q0)
% function recalq_two converts 4 market shares (including bundle goods) to
% 3 shares: og, two singletons

q01 = q0(1,1);
q02 = q0(2,1);
q03 = q0(3,1);
q04 = q0(4,1);

q1 = zeros(3,1);
q1(1,1) = q01;
q1(2,1) = q02 + q04;
q1(3,1) = q03 + q04;
end
