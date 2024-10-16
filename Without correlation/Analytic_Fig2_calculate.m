for i=1:length(nR)
    for Icase=1:3 %Quantization bit 지정, b에 따른 상수 지정, b infty인 경우 근사식 이용
        if Icase==1, q_gain=0.8825; beta=0.1175; 
        elseif Icase==2, q_gain=0.96546; beta=0.03454;
        else q_gain=0.990503; beta=0.009497;
        end

        G = @(z) exp(-z) .* (z.^nR(i)) .* (( gammainc(z, nR(i), 'lower')).^(nT - 1));
        G_q = @(z) exp(-z) .* (z.^(weibull_sum_parameter(1,i) + 1 / weibull_sum_parameter(2,i) - 1)) ...
            .* (( gammainc(z, weibull_sum_parameter(1,i), 'lower')).^(nT - 1));
        num_func=@(z) exp(-z) .* (z.^(nR(i)+1)) .* (( gammainc(z, nR(i), 'lower')).^(nT - 1));

        numerator=T_SNR*(1-beta)*quadgk(num_func, 0,  Inf, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        gaussian_noise=quadgk(@(z) G(z), 0,  Inf, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        quantization_noise=T_SNR*beta*quadgk(@(z) G_q(z), 0,  Inf, 'RelTol', 1e-20, 'AbsTol', 1e-20);
        if isnan(quantization_noise)
            denominator=gaussian_noise;%+2*T_SNR*beta*gamma(nR(i)+1);
        else
            denominator=gaussian_noise+quantization_noise;
        end

        R_analytic(Icase,i)=log2(1+(numerator/denominator));    
    end
end