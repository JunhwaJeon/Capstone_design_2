close all; clear; clc;

nT=32; nR=4;
v = [0, 0.2, 0.5, 0.8, 1]; % Teoplitz correlation coefficient 0.2
iter=1e3;
ind_err_vector=zeros(length(v), nT, nT);
ind_err_chi=zeros(length(v), nT, nT);

for Ccase=1:length(v)
    R_t = zeros(nT, nT);
    
    % R 행렬의 요소 계산
    for i = 1:nT
        for j = 1:nT
            R_t(i, j) = double(v(Ccase)^abs(i - j));
        end
    end
    phi_t=sqrtm(R_t); % Kroneker model antenna correlation matrix
    for i=1:iter
        W=sqrt(1/2)*complex(randn(nR, nT), randn(nR, nT));
        H=W*phi_t;
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
end
