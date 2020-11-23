function [ Sum_Rate, P, Theta ] = Opt_func( M,N,K,P_max,sigma,H,G )

%% init
Phi = diag(exp(1j*(rand(M, 1)*2*pi))); % phase shifters

%P = randn(N, K) + 1i*randn(N, K);
% ZF
H_all = H'*Phi*G;
P = H_all'*inv(H_all*H_all');
P = P/norm(P,'fro'); 

SumRate = [];
iter = 0;


%% Iterations UPDATE
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
    end
    Alpha(k) = S_gain/(IN_gain+sigma);
    sum_rate = sum_rate + log2(1+Alpha(k));
end

SumRate = [SumRate, sum_rate];



Beta = zeros(K,1);

for k=1:K
    S_gain1 = H_all(k,:)*P(:,k);
    IN_gain1 = 0;
    for j=1:K
        IN_temp1 = H_all(k,:)*P(:,j);
        IN_gain1 = IN_gain1 + abs(IN_temp1)^2;
    end
    Beta(k) = sqrt(1+Alpha(k))*S_gain1/(IN_gain1+sigma);
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


Theta = conj(diag(Phi));


V = cell(K,K);
for k=1:K
    v0 = diag(H(:,k)')*G;
    for j=1:K
        v = v0*P(:,j);
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


H_all = H'*Phi*G;



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

