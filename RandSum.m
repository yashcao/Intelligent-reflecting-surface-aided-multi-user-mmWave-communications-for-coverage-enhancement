function [ Sum_Rate ] = RandSum( bit, M, N,K,P_max,sigma,H,G  )
%RANDSUM 此处显示有关此函数的摘要
iter = 0;
Sum_rand = [];
%% init
Theta = exp(1j*(rand(M, 1)*2*pi));
if bit ~= 0
    Theta = Discerete(bit, Theta);
end

Phi = diag(Theta'); % phase shifters

% ZF
H_all = H'*Phi*G;
P = H_all'*inv(H_all*H_all');
P = P/norm(P,'fro'); 

while(1)
sum_rate = 0;
iter = iter+1;
% disp('iter:');fprintf('%c', 8);% 删掉换行符
% disp(iter);

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
        %IN_gain = IN_gain + IN_gain;
    end
    Alpha(k) = S_gain/(IN_gain+sigma);
    sum_rate = sum_rate + log2(1+Alpha(k));
end

Sum_rand = [Sum_rand, sum_rate];

Beta = zeros(K,1);
%f_max = 0;

for k=1:K
    S_gain1 = H_all(k,:)*P(:,k);
    IN_gain1 = 0;
    for j=1:K
        IN_temp1 = H_all(k,:)*P(:,j);
        IN_gain1 = IN_gain1 + abs(IN_temp1)^2;
        %IN_gain1 = IN_gain1 + real(h_12(k,:)*(P(:,j)*P(:,j)')*h_12(k,:)');
    end
    Beta(k) = sqrt(1+Alpha(k))*S_gain1/(IN_gain1+sigma);
    %f_max = f_max+2*sqrt(1+Alpha(k))*real(conj(Beta(k,:))*S_gain1)-Beta(k,:)'*Beta(k,:)*(IN_gain1+sigma);
end


mu = diag(sqrt(1+Alpha))*diag(Beta);
T = H_all'*diag(abs(Beta).^2)*H_all;

% syms  dual_v
power = @(dual_v) norm( inv(dual_v*eye(N)+T)*H_all'*mu, 'fro' )^2 - P_max;

low = 1e-6;
high = 50;
tolerance = 1e-5;
dual_v = Bisection(power, low, high, tolerance);

P = inv(dual_v*eye(N)+T)*H_all'*mu;


%% if converge
if iter >= 3
    diff = Sum_rand(iter) - Sum_rand(iter-1);
    if abs(diff)/Sum_rand(iter) <= 5e-3
        Sum_Rate = Sum_rand(end);
        break;
    end
end

end


end

