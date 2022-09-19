function [llobj,li] = loglikelihood_wog(theta,inp)

li = pyihat_wog(theta,inp);

logli = log(max(li,0.00000001));

llobj = -sum(logli);
end