function vars_matrix = mlevars(param0,inp)

stpsize = 0.00001;

np = inp.np;
n = inp.n;

deltalni_k = zeros(n,np);

pyihatmat = zeros(n,np);
pyihatmat0 = zeros(n,np);

for k =1:np
    param = param0;
    param(k,1) = param0(k,1)*(1+stpsize); % Change param values
    
    [~,pyihat0] = loglikelihood(param,inp);
    [~,pyihat00] = loglikelihood(param0,inp);
    pyihat = pyihat0';
    pyihat0 = pyihat00';
    
    pyihatmat(:,k) = pyihat;
    pyihatmat0(:,k) = pyihat0;
    
    deltalni = (pyihat - pyihat0)/(stpsize*param0(k,1));
    
    deltalni_k(:,k) = deltalni;
end

stack_indiv_var = zeros(np,np,n);

for i = 1:n
    stack_indiv_var(:,:,i) = deltalni_k(i,:)'*deltalni_k(i,:);
end

sum_indiv_var = sum(stack_indiv_var,3);

vars_matrix = inv(sum_indiv_var);
end


