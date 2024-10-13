close all; clear; clc;

R_t=[1, 0.2, 0.2, 0.2;...
        0.2, 1, 0.2, 0.2;...
        0.2, 0.2, 1, 0.2;...
        0.2, 0.2, 0.2, 1];
iter=10^6;
expectation_H=zeros(4,4);
for i=1:iter
    H=complex(randn(4,4),randn(4,4))*sqrt(1/2);
    H_corr=H*sqrt(R_t);
    expectation_H=expectation_H+H_corr;
end
expectation_H=expectation_H/iter;