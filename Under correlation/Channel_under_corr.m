close all; clear; clc;

R_t=[1, 0.2, 0.2, 0.2;...
        0.2, 1, 0.2, 0.2;...
        0.2, 0.2, 1, 0.2;...
        0.2, 0.2, 0.2, 1];
t_corr=sqrt(R_t);
iter=1e5;
expectation_H_corr=zeros(4,4);
expectation_H=zeros(4,4);
for i=1:iter
    H=complex(randn(4,4),randn(4,4))*sqrt(1/2);
    H_corr=H*t_corr;
    expectation_H=expectation_H+H;
    expectation_H_corr=expectation_H_corr+H_corr;
end
expectation_H_corr=expectation_H_corr/iter;
expectation_H=expectation_H/iter;