clear; clc;
syms z m k nt t

% 심볼릭 감마 함수 정의 (gammainc 대체)
assume(z > 0 & m > 0);
gamma_lower_sym = int(exp(-t) * t^(m-1), t, 0, z);

% g_q 함수 정의 (m, k, nt는 파라미터로 사용)
g_q = exp(-z) * z^(m + 1/k - 1) * (gamma_lower_sym^(nt - 1));

% z에 대한 부정적분 수행
F = int(g_q, z);

% 결과 출력
disp('g_q에 대한 부정적분 결과:');
disp(F);