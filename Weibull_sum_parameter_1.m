close all; clear; clc;
syms mu k

iid_weibull_moment=load('iid_weibull_moment.mat').iid_weibull_moment_124;
weibull_sum_mu_k=zeros(2,length(iid_weibull_moment(1,:)));

for i=1:length(iid_weibull_moment(1,:))
    eq1 = 2*gammaln(mu+1/k) - gammaln(mu+2/k) - gammaln(mu) == 2*log(iid_weibull_moment(1,i)) - log(iid_weibull_moment(2,i));
    eq2 = 2*gammaln(mu+2/k) - gammaln(mu+4/k) - gammaln(mu) == 2*log(iid_weibull_moment(2,i)) - log(iid_weibull_moment(3,i));

    % Provide initial guesses for mu and k
    initial_guess = [1, 1];  % Adjust the initial guess as necessary
    
    % Solve the system numerically with initial guesses
    [mu_sol, k_sol] = vpasolve([eq1, eq2], [mu, k], initial_guess);

    weibull_sum_mu_k(1,i)=mu_sol;
    weibull_sum_mu_k(2,i)=k_sol;
end

save('weibull_sum_parameter.mat','weibull_sum_mu_k'); % Nr=4,6,8, ... ,32 까지의 weibull sum parameter