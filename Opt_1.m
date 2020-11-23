function [ Sum_Rate, P, Theta ] = Opt_1( M,N,K,P_max,sigma,H,G )
%function [ SumRate, f4, Sum_dis, f5, iter, P, Theta ] = Opt_1( M,N,K,P_max,sigma,H,G )
%OPT_1 此处显示有关此函数的摘要
%% init

Phi = diag(exp(1j*(rand(M, 1)*2*pi))); % phase shifters

%P = randn(N, K) + 1i*randn(N, K);
% ZF
H_all = H'*Phi*G;
P = H_all'*inv(H_all*H_all');
P = P/norm(P,'fro'); 
% P = H_all';
% P = inv(H_all'*H_all+(sigma*K/P_max)*eye(N))*H_all';
% P = inv(H_all'*H_all)*H_all';

SumRate = [];
%loop = 0;
iter = 0;


%% Iterations UPDATE
%for iter=1:times
while(1)
iter = iter+1;
% disp('iter:');fprintf('%c', 8);% 删掉换行符
% disp(iter);


Alpha = zeros(K,1);
sum_rate = 0;


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

SumRate = [SumRate, sum_rate];



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
% P = H_all'*inv(H_all*H_all');
%P = P/norm(P,'fro'); 

Theta = conj(diag(Phi));
%Theta = reshape(Phi,[],1);


V = cell(K,K);
for k=1:K
    v0 = diag(H(:,k)')*G;
    for j=1:K
        v = v0*P(:,j);
        %v = reshape(v,[],1); % vec
        V(k,j) = {v};
    end
end


rho = zeros(K,1);

for k=1:K
    S_gain2 = sqrt(1+Alpha(k))*Theta'*V{k,k};
    IN_gain2 = 0;
    for j=1:K
        IN_temp2 = Theta'*V{k,j};
        IN_gain2 = IN_gain2 + abs(IN_temp2)^2;
        %IN_gain2 = IN_gain2 + real(Theta'*(V{k,j}*V{k,j}')*Theta);
    end
    rho(k) = S_gain2/(IN_gain2+sigma);
end

A = zeros(M, M);
b = zeros(M, 1);
%C = 0;
for k=1:K
    A = A + ([V{k,:}]*[V{k,:}]')*abs(rho(k))^2;
    b = b + sqrt(1+Alpha(k))*conj(rho(k))*V{k,k};
    %C = C + abs(rho(k))^2*sigma;
end

%-Theta'*A*Theta

% if(A==A') 
%     fprintf('是对称矩阵\n');
% else
%     fprintf('不是对称矩阵\n');
% end


cvx_begin sdp quiet
    variable f_SDP
    variable zeta(M)
    maximize (f_SDP-trace(diag(zeta)))
    subject to
        for m=1:M
            zeta(m) >= 0;
        end
        [A + diag(zeta)   b; ...
            b'   -f_SDP] >= 0;
cvx_end

% xxx = inv(A+diag(zeta))*(b*b')*inv(A+diag(zeta));
% diag(xxx)
Theta_opt = inv(A + diag(zeta))*b;

% 连续等幅度投影
Theta_opt = exp(1i*angle(Theta_opt));
% 离散等幅度投影
% bit = 1;
% Theta_opt1 = Discerete(bit, Theta_opt);
% bit = 2;
% Theta_opt2 = Discerete(bit, Theta_opt);

% abs(Theta_opt)
%f4(iter) = -real(Theta_opt'*A*Theta_opt) + 2*real(Theta_opt'*b) - C;


Phi = diag(Theta_opt');
% Phi = diag(Theta);



H_all = H'*Phi*G;
% P = H_all'*inv(H_all*H_all');
% P = H_all';
% P = P/norm(P,'fro'); 


% if converge
% if iter >= 3
%     diff = SumRate(iter) - SumRate(iter-1);
%     if abs(diff) <= 0.01
%         loop = loop+1;
%         
%         if loop >= 5
%             Sum_Rate = SumRate(end);
%             break;
%         end
%     else
%         loop = 0;
%     end
% end

%% if converge
if iter >= 3
    diff = SumRate(iter) - SumRate(iter-1);
    if abs(diff)/SumRate(iter) <= 1e-2
        Sum_Rate = SumRate(end);
        break;
    end
end    
    

end


end

