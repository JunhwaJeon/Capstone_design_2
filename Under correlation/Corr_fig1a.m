clear; clc; close all;
%% 4*4 MIMO에서 우선 실험

%% 파라미터 설정
T_snr_db=-15:1:20;
T_snr=10.^(T_snr_db./10);
iter=1e3;
R_t=[1, 0.2, 0.2, 0.2;... % transmit antenna correlation matrix
        0.2, 1, 0.2, 0.2;...
        0.2, 0.2, 1, 0.2;...
        0.2, 0.2, 0.2, 1]; 
phi_t=sqrt(R_t); % Kroneker model antenna correlation matrix
nT=4; nR=4;
R=zeros(5, length(T_snr));
R_randi=zeros(nT);

for Qcase=1:5
    if Qcase==1, q_gain=0.6364; beta=0.3636; 
    elseif Qcase==2, q_gain=0.8825; beta=0.1175; 
    elseif Qcase==3, q_gain=0.96546; beta=0.03454;
    elseif Qcase==4, q_gain=0.990503; beta=0.009497;
    else q_gain=1; beta=0; 
    end
    for p=1:length(T_snr)
        for i=1:iter
            H=sqrt(1/2)*complex(randn(4,4), randn(4,4));
            for j=1:nT % Transmit Antenna Selection
                sum_four_sqr=0;
                norm_sqr=(H(:,j))'*H(:,j);
                for k=1:nR
                    sum_four_sqr=sum_four_sqr+abs(H(k,j))^4;
                end
                R_candi(j)=log2(1+(T_snr(p)*q_gain*(norm_sqr)^2)/(norm_sqr+T_snr(p)*beta*sum_four_sqr));
            end
            R(Qcase, p)=R(Qcase, p)+max(R_candi);
        end
    end
end

R=R/iter;
plot(T_snr_db,R(1,:),'b-', T_snr_db,R(2,:),'b-', T_snr_db,R(3,:),'b-');
hold on, grid on,
plot(T_snr_db,R(4,:),'b-', T_snr_db, R(5,:),'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-15,21]);
legend('Simulation rate in Eq.(6)')