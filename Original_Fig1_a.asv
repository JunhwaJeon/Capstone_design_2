%원 논문 Fig1.(a) 그려보기
clear all, close all
clc;

%%파라미터 설정
T_SNR_dB=-15:1:20; %SNR 범위 설정
T_SNR_linear=10.^(T_SNR_dB/10); %linear 스케일 SNR설정
N_iter=1000; %반복 횟수 (Ergodic capacity 구하기 위해서) 
sq2 = sqrt(0.5); %상수 지정
nT=32; nR=8; %MIMO Scale 지정
analytic=load("Fig_1a_analytic.mat").R;

%%Quantization bit 지정, b에 따른 상수 지정, b infty인 경우 근사식 이용
for Icase=1:5
    if Icase==1, q_gain=0.6364; beta=0.3634; 
    elseif Icase==2, q_gain=0.8825; beta=0.1175; 
    elseif Icase==3, q_gain=0.96546; beta=0.03454;
    elseif Icase==4, q_gain=0.990503; beta=0.009497;
    else q_gain=1; beta=0; 
    end

R(Icase,:) = zeros(1,length(T_SNR_dB));%Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)
R_candi=linspace(0,0,nT); %Maximum 선택 위한 후보값 담을 행렬(벡터) 지정

%% Ergodic Capacity 계산
for i=1:length(T_SNR_dB)
    for iter=1:N_iter %반복
        H = sq2*complex(randn(nR,nT),randn(nR,nT)); %Complex Circular Gaussian channel (Rayleigh)
        sum_four_sqr=0; norm_sqr=0;
        for j=1:nT %Transmit Antenna Selection
            norm_sqr=(H(:,j))'*H(:,j);
            sum_four_sqr=0;
            for k=1:nR
                sum_four_sqr=sum_four_sqr+abs(H(k,j))^4;
            end
            R_candi(j)=log2(1+(T_SNR_linear(i)*(1-beta)*(norm_sqr)^2)/(norm_sqr+T_SNR_linear(i)*(beta)*sum_four_sqr));
        end
        R(Icase,i)=R(Icase,i)+max(R_candi);
    end
end
end

R = R/N_iter; %Expectation 계산

plot(T_SNR_dB,R(1,:),'b-', T_SNR_dB,R(2,:),'b-', T_SNR_dB,R(3,:),'b-', T_SNR_dB,R(4,:),'b-', T_SNR_dB, R(5,:),'b-');
hold on, grid on,
plot(T_SNR_dB, analytic(1,:),'ro', T_SNR_dB, analytic(2,:),'ro', T_SNR_dB, analytic(3,:),'ro', T_SNR_dB, analytic(4,:),'ro', T_SNR_dB, analytic(5,:),'ro');
xlabel('Transmit SNR[dB]'); ylabel('Ergodic Rate [bps/Hz]');