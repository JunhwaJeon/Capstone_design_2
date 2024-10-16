close all; clear; clc;
syms mu k

iid_weibull_moment=load('iid_weibull_moment_rev1.mat').iid_weibull_moment_124;
weibull_sum_mu_k=zeros(2,length(iid_weibull_moment(1,:)));

% 파라미터 범위 설정 (예: 양수)
lb = [0, 0];  % 하한값
ub = [Inf, Inf];  % 상한값

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-6);

for i=1:length(iid_weibull_moment(1,:))
    M1 = iid_weibull_moment(1, i);
    M2 = iid_weibull_moment(2, i);
    M3 = iid_weibull_moment(3, i);
    
    % 목적 함수 정의
    fun = @(x) [
        2*gammaln(x(1) + 1 / x(2)) - gammaln(x(1) + 2 / x(2)) - gammaln(x(1)) - log(M1^2 / M2);
        2*gammaln(x(1) + 2 / x(2)) - gammaln(x(1) + 4 / x(2)) - gammaln(x(1)) - log(M2^2 / M3);
    ];

    % 초기 추정값 설정 (양수 값)
    x0 = [5, 1];

    % 비선형 제약 조건 정의: gammaln의 입력이 양수가 되도록 함
    nonlcon = @(x) deal([], [-(x(1) + 1/x(2)); -(x(1) + 2/x(2)); -(x(1) + 4/x(2))]);

    % 수치적으로 방정식 풀이
    [x_sol, fval, exitflag] = fmincon(@(x) sum(fun(x).^2), x0, [], [], [], [], lb, ub, nonlcon, options);

    if exitflag > 0 && x_sol(1) > 0 && x_sol(2) > 0
        weibull_sum_mu_k(1, i) = x_sol(1);
        weibull_sum_mu_k(2, i) = x_sol(2);
    else
        weibull_sum_mu_k(1, i) = NaN;
        weibull_sum_mu_k(2, i) = NaN;
        warning('Iteration %d: Solution not found or invalid.', i);
    end
end

save('weibull_sum_parameter_rev1.mat','weibull_sum_mu_k'); % Nr=4,6,8, ... ,32 까지의 weibull sum parameter