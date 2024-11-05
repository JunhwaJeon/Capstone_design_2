clear; clc; close all;

%% 파라미터 설정
T_snr_db=-20:5:30;
T_snr=10.^(T_snr_db./10);
iter=1e3;
nT=32; nR=1;
v = [0, 0.2, 0.5, 0.8, 1]; % Teoplitz correlation coefficient 0.2
R_corr=zeros(length(v), 5, length(T_snr));
R_anal = zeros(5,length(T_snr)); % Analytic Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)

for Ccase=1:length(v)
    R_t = zeros(nT);
    B=zeros(1,nT);

    % R 행렬의 요소 계산
    for i = 1:nT
        for j = 1:nT
            R_t(i, j) = double(v(Ccase)^abs(i - j));
        end
    end
    phi_t=sqrtm(double(R_t)); % Kroneker model antenna correlation matrix
    for i = 1:nT
        B(i)=sum(phi_t(:,i).^2);
    end

for Qcase=1:5
    if Qcase==1, q_gain=0.6364; beta=0.3636; 
    elseif Qcase==2, q_gain=0.8825; beta=0.1175; 
    elseif Qcase==3, q_gain=0.96546; beta=0.03454;
    elseif Qcase==4, q_gain=0.990503; beta=0.009497;
    else q_gain=1; beta=0; 
    end

    for p=1:length(T_snr)
        for i=1:iter
            H_w=sqrt(1/2)*complex(randn(nR, nT), randn(nR, nT));
            H=H_w*phi_t;
            R_randi=zeros(nT);
            for j=1:nT % Transmit Antenna Selection
                sum_four_sqr=0;
                norm_sqr=(H(:,j))'*H(:,j);
                norm_sqr_nocorr=(H_w(:,j))'*H_w(:,j);
                norm_sqr_cal=0;
                for k=1:nR
                    norm_sqr_cal=norm_sqr_cal+phi_t(:,j)'*H_w(k,:)'*H_w(k,:)*phi_t(:,j);
                end

                for k=1:nR
                    sum_four_sqr=sum_four_sqr+abs(H(k,j))^4;
                end
                R_candi(j)=log2(1+(T_snr(p)*q_gain*(norm_sqr)^2)/(norm_sqr+T_snr(p)*beta*sum_four_sqr));
            end
            R_corr(Ccase, Qcase, p)=R_corr(Ccase, Qcase, p)+max(R_candi);
        end
    end

    for p=1:length(T_snr)
        
    end
end
end

R_corr=R_corr/iter;

figure('Name', 'Correlation coefficient=0.2', 'Numbertitle', 'off');
plot(T_snr_db, reshape(R_corr(1,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(1,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(1,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(1,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(1,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,20]);
legend('Simulation rate in Eq.(6)')

figure('Name', 'Correlation coefficient=0.5', 'Numbertitle', 'off' )
plot(T_snr_db, reshape(R_corr(2,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(2,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(2,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(2,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(2,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,20]);
legend('Simulation rate in Eq.(6)')

figure('Name', 'Correlation coefficient=0.8', 'Numbertitle', 'off')
plot(T_snr_db, reshape(R_corr(3,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(3,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(3,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(3,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(3,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,20]);
legend('Simulation rate in Eq.(6)')