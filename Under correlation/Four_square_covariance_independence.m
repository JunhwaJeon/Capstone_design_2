close all; clear; clc;

nT=32; nR=4;
v = [0, 0.2, 0.5, 0.8, 1]; % Teoplitz correlation coefficient
iter=5*1e3;
ind_vector=zeros(length(v), nT, nT);
ind_err_vec=zeros(length(v), nT, nT);
Real_val_jj=zeros(length(v), nT, nT);

for Ccase=1:length(v)
    R_t = zeros(nT, nT);
    
    %% Making Kroneker model transmit antenna correlation matrix
    for i = 1:nT
        for j = 1:nT
            R_t(i, j) = double(v(Ccase)^abs(i - j));
        end
    end
    phi_t=sqrtm(R_t); % Kroneker model antenna correlation matrix

    %% Compute Joint expectation of Sum of four squares
    
    %{
    for i=1:iter
        W=sqrt(1/2)*complex(randn(nR, nT), randn(nR, nT)); % uncorrelated channel
        H=W*phi_t; % correlated channel
        for jj1=1:nT
            for jj2=1:nT
                if_independent_vector=0;
                if_independent_chi=0;
                Z_jj1=0; Z_jj2=0;
                for k=1:nR
                    Z_jj1=Z_jj1+abs(H(k, jj1))^4; % Real Value of the sum of transmit antenna
                    Z_jj2=Z_jj2+abs(H(k, jj2))^4;
                end
                non_independent=Z_jj1*Z_jj2;
                for kk1=1:nR
                    for kk2=1:nR
                        if_independent_chi=if_independent_chi+(chi2rnd(2)*chi2rnd(2))*...
                            phi_t(:,jj1)'*W(kk1,:)'*W(kk1,:)*phi_t(:, jj1)*phi_t(:,jj2)'*W(kk2,:)'*W(kk2,:)*phi_t(:,jj2);
                        if_independent_vector=if_independent_vector+W(kk1,:)*phi_t(:, jj1)*phi_t(:,jj1)'*W(kk1,:)'*...
                            W(kk2,:)*phi_t(:,jj2)*phi_t(:,jj2)'*W(kk2,:)'*...
                        phi_t(:,jj1)'*W(kk1,:)'*W(kk1,:)*phi_t(:, jj1)*phi_t(:,jj2)'*W(kk2,:)'*W(kk2,:)*phi_t(:,jj2);
                    end
                end
            ind_err_vector(Ccase, jj1, jj2)=abs(non_independent-if_independent_vector)/non_independent;
            ind_err_chi(Ccase, jj1, jj2)=abs(non_independent-if_independent_chi)/non_independent;
            end
        end
    end
    ind_err_vector=ind_err_vector/iter;
    ind_err_chi=ind_err_chi/iter;
    %}
    for jj1=1:nT
        for jj2=1:nT
            A_kk=zeros(nR, nR); Square_jnt_exp=zeros(nR, nR);
            Z_jj1=zeros(nR, nR); Z_jj2=zeros(nR, nR);
            Real_val=zeros(nR, nR);
            for i=1:iter
                W=sqrt(1/2)*complex(randn(nR, nT), randn(nR, nT)); % uncorrelated channel
                H=W*phi_t; % correlated channel
                for kk1=1:nR
                    for kk2=1:nR
                        A_kk(kk1, kk2)=A_kk(kk1, kk2)+...
                            chi2rnd(2)*chi2rnd(2);
                            %W(kk1,:)*phi_t(:, jj1)*phi_t(:,jj1)'*W(kk1,:)'*W(kk2,:)*phi_t(:,jj2)*phi_t(:,jj2)'*W(kk2,:)';
                        Square_jnt_exp(kk1, kk2)=Square_jnt_exp(kk1, kk2)+abs(W(kk1,:)*phi_t(:, jj1))^2*abs(W(kk2,:)*phi_t(:,jj2))^2;
                        Real_val(kk1, kk2)=Real_val(kk1, kk2)+(abs(H(kk1, jj1))^4)*(abs(H(kk2, jj2))^4);
                    end
                end
            end
            A_kk=A_kk/iter; Square_jnt_exp=Square_jnt_exp/iter; Real_val=Real_val/iter;
            If_ind_mtx=A_kk.*Square_jnt_exp;
            If_ind=sum(If_ind_mtx, 'all');
            Non_ind=sum(Real_val, 'all');
            Real_val_jj(Ccase, jj1, jj2)=Non_ind;
            ind_vector(Ccase, jj1, jj2)=If_ind;
        end
    end
end
ind_err_vec=abs(Real_val_jj-ind_vector)./Real_val_jj;
off_diag_err=zeros(length(v));
for i=1:length(v)
    off_diag_err(i)=(sum(squeeze(ind_err_vec(i,:,:)), 'all')-sum(diag(squeeze(ind_err_vec(i,:,:))),'all'))/(32*31);
end


