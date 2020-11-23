function [ DisSumRate ] = DiscreteSum( bit,Theta_opt,K,P,sigma,H,G )

sum_rate = 0;

Theta_opt = Discerete(bit, Theta_opt);
Phi = diag(Theta_opt');

H_all = H'*Phi*G;

for k=1:K
    S_temp = H_all(k,:)*P(:,k);
    S_gain = abs(S_temp)^2;
    IN_gain = 0;

    for j=1:K
        if j==k
            continue;
        end
        IN_temp = H_all(k,:)*P(:,j);
        IN_gain = IN_gain + abs(IN_temp)^2;
    end
    Alpha(k) = S_gain/(IN_gain+sigma);
    sum_rate = sum_rate + log2(1+Alpha(k));
end

DisSumRate = sum_rate;

end

