clear; clc;
syms z m1 m2 k nt t

% 심볼릭 감마 함수 정의 (gammainc 대체)
assume(z > 0 & m1 > 0 & m2>0);
gamma_upper_sym1 = int(exp(-t) * t^(m1-1), t, z, Inf);
gamma_lower_sym2 = int(exp(-t) * t^(m2-1), t, 0, z);

% g_q 함수 정의 (m, k, nt는 파라미터로 사용)
g_q = gamma_upper_sym1.* gamma_lower_sym2;

% z에 대한 부정적분 수행
F = int(g_q, z, 0, Inf);

% 결과 출력
disp('g_q에 대한 부정적분 결과:');
disp(F);