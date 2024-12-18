function iid_weibull_moment=Weibull_parameter(Nr_array, Bi ) % Nr_array: Nr 개수를 가지고 있는 배열, Bi: Scaling factor, Tx에 따라 달라짐
iid_weibull_moment=zeros(3,length(Nr_array));
    for m=1:length(Nr_array)
        iid_weibull_moment(1,m) = computeExpectationR(1, Nr_array(m), Bi);
        iid_weibull_moment(2,m) = computeExpectationR(2, Nr_array(m), Bi);
        iid_weibull_moment(3,m) = computeExpectationR(4, Nr_array(m), Bi);
    end

    %% Parameters only dependes on the number of receive antennas
    function E_R_n = computeExpectationR(n, Nr, Bi)
    
        E_R_n=0;
        all_combinations=generate_sequences(n, Nr);
    
        for j=1:length(all_combinations(:,1))
            all_combinations(j,:)=sort(all_combinations(j,:),"descend");
        end
        for j=1:length(all_combinations(:,1))
            vector=all_combinations(j,:);
            part=1;
            for k=1:length(vector)-1
            part=part*nchoosek(vector(k),vector(k+1))*gamma(1+2*(vector(k)-vector(k+1))/Bi);
            end
            E_R_n=E_R_n+part*gamma(1+2*vector(length(vector))/Bi);
        end
        
        function result = generate_sequences(n, Nr)
        % Nr 길이의 벡터 생성, 첫 번째 요소는 항상 n, 이후 값은 같거나 감소
        result = generate_recursive(n, Nr-1); % 첫 번째 요소 이후 비내림차순 시퀀스 생성
        result = [n * ones(size(result, 1), 1), result]; % 첫 번째 요소를 n으로 설정
    end
    
    function sequences = generate_recursive(max_value, length)
        % 재귀적으로 같거나 감소하는 벡터 생성
        if length == 0
            sequences = [];
        elseif length == 1
            sequences = (0:max_value)'; % 마지막 요소는 0에서 max_value까지
        else
            sequences = [];
            for i = max_value:-1:0
                sub_sequences = generate_recursive(i, length-1);
                sequences = [sequences; i * ones(size(sub_sequences, 1), 1), sub_sequences];
            end
        end
    end

end
end

