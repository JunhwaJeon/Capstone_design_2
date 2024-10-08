close all; clear; clc;

%% 실험환경
% under fixed number of transmit and receive antennas at 32, 8
% Calculate analatical value of Transmit SNR vs. Ergodic rate

%% 파라미터 설정
T_SNR_dB=-15:1:20; %SNR 범위 설정
T_SNR=10.^(T_SNR_dB/10); %linear 스케일 SNR설정
nT=32; nR=8; %MIMO Scale 지정
q_gain=0;
R = zeros(5,length(T_SNR_dB)); % Analatic Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)
weibull_sum_parameter=load('weibull_sum_parameter.mat').weibull_sum_mu_k;
weibull_sum_parameter=weibull_sum_parameter(:,3);
G=@(z) exp(-z).*(z.^nR).*(gammainc(z,nR).^(nT-1));
G_q=@(z) (2*gamma(nR+1)/gamma(weibull_sum_parameter(1)+1/weibull_sum_parameter(2))).*exp(-z).* ...
    (z.^(weibull_sum_parameter(1)+1/weibull_sum_parameter(2)-1)).*(gammainc(z, weibull_sum_parameter(1)).^(nT-1));

%% Quantization bit 지정, b에 따른 상수 지정, b infty인 경우 근사식 이용
for Qcase=1:5
    if Qcase==1, q_gain=0.6364; beta=0.3636; 
    elseif Qcase==2, q_gain=0.8825; beta=0.1175; 
    elseif Qcase==3, q_gain=0.96546; beta=0.03454;
    elseif Qcase==4, q_gain=0.990503; beta=0.009497;
    else q_gain=1; beta=0; 
    end
  
    for p=1:length(T_SNR)
        num_func=@(z) exp(-z).*(z.^(nR+1)).*(gammainc(z,nR).^(nT-1));
        numerator=T_SNR(p)*(1-beta)*integral(num_func, 0, Inf);
        denominator=integral(G, 0, Inf)+T_SNR(p)*beta*integral(G_q, 0, Inf);
        R(Qcase,p)=log2(1+numerator/denominator);
    end
end
save('Fig_1a_analatic.mat','R');

plot(T_SNR_dB,R(1,:),'ro', T_SNR_dB,R(2,:),'ro', T_SNR_dB,R(3,:),'ro');
hold on, grid on,
plot(T_SNR_dB,R(4,:),'ro', T_SNR_dB, R(5,:),'ro');
xlabel('Transmit SNR[dB]'); ylabel('Ergodic Rate [bps/Hz]');