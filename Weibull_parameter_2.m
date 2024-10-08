close all; clear; clc;

iid_weibull_moment_124=zeros(3,15);

for i=1:15
    iid_weibull_moment_124(1,i)=computeExpectationR(1, 2*i+2);
    iid_weibull_moment_124(2,i)=computeExpectationR(2, 2*i+2);
    iid_weibull_moment_124(3,i)=computeExpectationR(4, 2*i+2);
end

save('iid_weibull_moment.mat',"iid_weibull_moment_124");

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