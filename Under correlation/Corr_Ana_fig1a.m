clear; clc; close all;
syms z
%% 파라미터 설정
T_snr_db=-15:1:20;
T_snr=10.^(T_snr_db./10);
nT=32; nR=8;
v = [0.2, 0.4, 0.6]; % Teoplitz correlation coefficient 0.2
R_anal = zeros(length(v),5,length(T_snr)); % Analatic Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)

for Ccase=1:length(v) % Correlation Case 0.2, 0.5, 0.8
    R_t = zeros(nT);
    B=zeros(1,nT);
    % R 행렬의 요소 계산
    for i = 1:nT
        for j = 1:nT
            R_t(i, j) = v(Ccase)^abs(i - j);
        end
    end
    phi_t=sqrtm(R_t); % Kroneker model antenna correlation matrix
    for i = 1:nT
        B(i)=sum(phi_t(:,i).^2);
    end

    % G(z)와 G_q(z)를 수치적으로 정의하기
    G = @(z) 1; % G 함수 초기화
    for i = 1:nT
        G_i = @(z) gammainc(z / (2 * B(i)), nR, 'lower') / factorial(nR - 1);
        G = @(z) G(z) .* G_i(z); % 각 B(i)에 대해 곱셈 수행
    end

    % G_q(z) 정의
    G_q = @(z) 1; % G_q 함수 초기화
    for i = 1:nT
        sum_parameter = Weibull_sum_parameter(nR, B(i));
        moment = Weibull_parameter(nR, B(i));
        G_q_i = @(z) gammainc(((gamma(1 / sum_parameter(2)) .* z) / (beta(sum_parameter(1), 1 / sum_parameter(2)) * moment(1))).^sum_parameter(2), sum_parameter(1), 'lower');
        G_q = @(z) G_q(z) .* G_q_i(z); % 각 B(i)에 대해 곱셈 수행
    end

    % 적분을 위한 함수 정의 및 수치 적분 수행
    integrand_sig = @(z) (z.^2) .* G(z);
    integrand_noise = @(z) z .* G(z);
    integrand_qant_noise = @(z) z .* G_q(z);

    % 수치적 적분 수행
    Sig_term = integral(integrand_sig, 0, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-12);
    Noise_term = integral(integrand_noise, 0, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-12);
    Qant_noise = integral(integrand_qant_noise, 0, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-12);


    for Qcase=1:5
        if Qcase==1, q_gain=0.6364; 
        elseif Qcase==2, q_gain=0.8825; 
        elseif Qcase==3, q_gain=0.96546; 
        elseif Qcase==4, q_gain=0.990503; 
        else q_gain=1; 
        end
        for p=1:length(T_snr)
            R_anal(Ccase, Qcase, p)=log2(1+(T_snr(p)*q_gain*Sig_term)/(Noise_term+T_snr(p)*(1-q_gain)*Qant_noise));
        end
    end
end
