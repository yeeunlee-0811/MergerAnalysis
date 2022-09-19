function [llobj,li] = loglikelihood(theta,inp)

li = pyihat(theta,inp);

logli = log(max(li,0.00000001));

llobj = -sum(logli);
end