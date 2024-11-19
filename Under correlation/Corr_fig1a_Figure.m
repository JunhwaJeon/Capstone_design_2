T_snr_db=-20:5:30;


figure('Name', 'Correlation coefficient=0', 'Numbertitle', 'off');
plot(T_snr_db, reshape(R_corr(1,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(1,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(1,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(1,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(1,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,30]);
legend('Simulation rate in Eq.(6)')

figure('Name', 'Correlation coefficient=0.2', 'Numbertitle', 'off');
plot(T_snr_db, reshape(R_corr(2,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(2,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(2,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(2,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(2,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,30]);
legend('Simulation rate in Eq.(6)')

figure('Name', 'Correlation coefficient=0.5', 'Numbertitle', 'off' )
plot(T_snr_db, reshape(R_corr(3,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(3,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(3,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(3,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(3,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,30]);
legend('Simulation rate in Eq.(6)')

figure('Name', 'Correlation coefficient=0.8', 'Numbertitle', 'off')
plot(T_snr_db, reshape(R_corr(4,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(4,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(4,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(4,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(4,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,20]);
legend('Simulation rate in Eq.(6)')

figure('Name', 'Correlation coefficient=1', 'Numbertitle', 'off')
plot(T_snr_db, reshape(R_corr(5,1,:),[],1) ,'b-', T_snr_db,reshape(R_corr(5,2,:),[],1) ,'b-', T_snr_db,reshape(R_corr(5,3,:),[],1) ,'b-');
hold on, grid on,
plot(T_snr_db,reshape(R_corr(5,4,:),[],1) ,'b-', T_snr_db, reshape(R_corr(5,5,:),[],1) ,'b-');
xlabel('Transmit SNR [dB]'); ylabel('Ergodic Rate [bps/Hz]');
xlim([-20,30]);
legend('Simulation rate in Eq.(6)')