clear all; close all; clc;
syms mu k

% Define the system of equations
eq1 = 2*gammaln(mu+1/k) - gammaln(mu+2/k) - gammaln(mu) == 2*log(computeExpectationR(1, 16)) - log(computeExpectationR(2, 16));
eq2 = 2*gammaln(mu+2/k) - gammaln(mu+4/k) - gammaln(mu) == 2*log(computeExpectationR(2, 16)) - log(computeExpectationR(4, 16));

% Provide initial guesses for mu and k
initial_guess = [1, 1];  % Adjust the initial guess as necessary

% Solve the system numerically with initial guesses
[mu_sol, k_sol] = vpasolve([eq1, eq2], [mu, k], initial_guess);

% Display the solutions
disp('Solutions for mu:');
disp(mu_sol);
disp('Solutions for k:');
disp(k_sol);

%% Parameters only dependes on the number of receive antennas
function E_R_n = computeExpectationR(n, M)
    % n: 차수
    % M: 차원의 개수
    E_R_vals = repmat({@(n) gamma(1 + 2*n)}, 1, M);

    E_R_n = 0;  % 결과를 저장할 변수
    
    % 재귀적으로 M-1 차원까지 합을 계산하는 함수
    function recursiveSum(k, indices, sum_term)
        if k == M-1
            % 마지막 차원일 경우
            n_remaining = n - sum(indices);  % 남은 n 값 계산
            if isempty(indices)  % 첫 번째 차원에서는 indices가 비어 있음
                binom_coeff = nchoosek(n, n_remaining);  % nchoosek(N, K) 호출 전에 값 확인
            else
                binom_coeff = 1;  % 이항 계수 초기값
                for j = 1:length(indices)
                    N = n - sum(indices(1:j-1));  % 남은 N 값 계산
                    K = indices(j);  % K 값
                    if K <= N && K >= 0  % nchoosek가 유효한지 확인
                        binom_coeff = binom_coeff * nchoosek(N, K);
                    else
                        binom_coeff = 0;  % 유효하지 않은 경우 이항계수는 0
                        break;
                    end
                end
            end
            % 남은 값과 함께 최종 계산
            if n_remaining >= 0
                E_R_n = E_R_n + sum_term * binom_coeff * E_R_vals{M}(n_remaining);
            end
        else
            % 다음 차원으로 이동
            for j = 0:n-sum(indices)
                recursiveSum(k+1, [indices, j], sum_term * E_R_vals{k}(n - sum(indices)));
            end
        end
    end

    % 첫 번째 차원부터 재귀 합을 시작
    recursiveSum(1, [], 1);
end