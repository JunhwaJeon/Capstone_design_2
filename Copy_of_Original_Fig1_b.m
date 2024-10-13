clear all, close all
clc;

%%파라미터 설정
T_SNR_dB=10; %SNR 범위 설정
T_SNR_linear=10.^(T_SNR_dB/10); %linear 스케일 SNR설정Line 6 T_SNR_linear=10.^(T_SNR_dB/10); %linear 스케일 SNR설정
N_iter=1000; %반복 횟수 (Ergodic capacity 구하기 위해서) 
sq2 = sqrt(0.5); %상수 지정
nT=32; 
nR=[4:1:32]; %MIMO Scale 지정
q_gain=0;

R_candi=linspace(0,0,nT); %Maximum 선택 위한 후보값 담을 행렬(벡터) 지정
R(:,:) = zeros(3,length(nR));%Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)
%% Ergodic Capacity 계산
for i=1:length(nR)
    for Icase=1:3 %Quantization bit 지정, b에 따른 상수 지정, b infty인 경우 근사식 이용
        if Icase==1, q_gain=0.8825; 
        elseif Icase==2, q_gain=0.96546;
        else q_gain=0.990503; 
        end
        for iter=1:N_iter %반복
            H = sq2*(randn(nR(i),nT)+1j*randn(nR(i),nT)); %Complex Circular Gaussian channel (Rayleigh)
            if nR(i)>=nT, HH = H'*H; else HH = H*H'; end
            for j=1:nT %Transmit Antenna Selection
            sum_four_sqr=0; norm_sqr=0;
            norm_sqr=(H(:,j))'*H(:,j);
                for k=1:nR(i)
                sum_four_sqr=sum_four_sqr+abs(H(k,j))^4;
                end
            R_candi(j)=log2(1+(T_SNR_linear*q_gain*(norm_sqr)^2)/(norm_sqr+T_SNR_linear*(1-q_gain)*sum_four_sqr));
            end
            R(Icase,i)=R(Icase,i)+max(R_candi);
        end
    end
end
R = R/N_iter; %Expectation 계산

plot(nR,R(1,:),'b-', nR,R(2,:),'b-', nR,R(3,:),'b-');
hold on, grid on,
xlabel('Number of Receive Antennas Nr'); ylabel('Ergodic Rate [bps/Hz]');
