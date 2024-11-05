
function sum_parameter=Weibull_sum_parameter(Nr_array, Bi)
    iid_parameter=Weibull_parameter(Nr_array, Bi);
    sum_parameter=zeros(2, length(Nr_array));
    syms mu k
    % 양수 조건 설정
    assume(mu > 0);
    assume(k > 0);

    for i = 1:length(Nr_array)
        M1 = iid_parameter(1, i);
        M2 = iid_parameter(2, i);
        M4 = iid_parameter(3, i);

        eq1 = 2*gammaln(mu + 1/k) - gammaln(mu + 2/k) - gammaln(mu) == 2*log(M1) - log(M2);
        eq2 = 2*gammaln(mu + 2/k) - gammaln(mu + 4/k) - gammaln(mu) == 2*log(M2) - log(M4);

        % Provide initial guesses for mu and k
        initial_guess = [1, 1];  % Adjust the initial guess as necessary

        % Solve the system numerically with initial guesses
        [mu_sol, k_sol] = vpasolve([eq1, eq2], [mu, k], initial_guess);

        % 만약 양수 해가 아닐 경우 처리
        if isempty(mu_sol) || isempty(k_sol) || mu_sol <= 0 || k_sol <= 0
            warning('No valid positive solution found for iteration %d', i);
            sum_parameter(1, i) = NaN; % 양수 해가 없는 경우 NaN으로 설정
            sum_parameter(2, i) = NaN;
        else
            sum_parameter(1, i) = double(mu_sol);
            sum_parameter(2, i) = double(k_sol);
        end
    end
end