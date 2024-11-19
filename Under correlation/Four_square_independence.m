close all; clear; clc;

nT=32; nR=4;
v = [0, 0.2, 0.5, 0.8, 1]; % Teoplitz correlation coefficient 0.2
iter=1e3;
ind_err_vector=zeros(length(v), nT);
ind_err_chi=zeros(length(v), nT);

for Ccase=1:length(v)
    R_t = zeros(nT, nT);
    
    % R 행렬의 요소 계산
    for i = 1:nT
        for j = 1:nT
            R_t(i, j) = v(Ccase)^abs(i - j);
        end
    end
    phi_t=sqrtm(R_t); % Kroneker model antenna correlation matrix
    for i=1:iter
        W=sqrt(1/2)*complex(randn(nR, nT), randn(nR, nT));
        H=W*phi_t;
        for j=1:nT
            non_independent=0;
            if_independent_vector=0;
            if_independent_chi=0;
            for k=1:nR
                non_independent=non_independent+abs(H(k,j))^4;
                if_independent_vector=if_independent_vector+ W(k,:)*phi_t(:,j)*phi_t(:,j)'*W(k,:)'*...
                    phi_t(:,j)'*W(k,:)'*W(k,:)*phi_t(:,j);
                if_independent_chi=if_independent_chi+chi2rnd(2)*phi_t(:,j)'*W(k,:)'*W(k,:)*phi_t(:,j);
            end
            ind_err_vector(Ccase, j)=abs(non_independent-if_independent_vector)/non_independent;
            ind_err_chi(Ccase, j)=abs(non_independent-if_independent_chi)/non_independent;
        end
    end
    ind_err_vector=ind_err_vector/iter;
    ind_err_chi=ind_err_chi/iter;
end
