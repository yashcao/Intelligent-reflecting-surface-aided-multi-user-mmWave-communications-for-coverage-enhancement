function [H2] = IRS_user(Mx, My, K)
% IRS-user channel matrix
M  = Mx*My;
H2 = zeros(M, K);


radis = 8; % 8m


Range = rand(K, 1) * (radis-1)+1;
a1 = (rand(K, 1)-0.5)*pi;
a2 = (rand(K, 1)-0.5)*pi;
user_pos = zeros(K, 3);

user_pos(:, 3) = Range.* cos(a1);
user_pos(:, 1) = Range.* (cos(a2).*sin(a1));
user_pos(:, 2) = Range.* (sin(a2).*sin(a1));



%% array steering vectors
%dd = 0.5;        % dd = d/lambda; d=lamda/2
lamda = 0.0107;  % 28 GHZ
d = lamda/2;

dx = (-Mx/2:Mx/2-1)*d+0.5*d;
dy = (-My/2:My/2-1)*d+0.5*d;


for k=1:K
    dis_0(k) = sqrt(sum([0, 0, 0] - user_pos(k)).^2);
end


for k=1:K
    Arr_U = zeros(Mx, My);
    
    for mx = 1:Mx
        for my = 1:My
            dis_P(mx, my) = sqrt(sum([dx(mx), dy(my), 0] - user_pos(k)).^2);
            Arr_U(mx, my) = (dis_0(k)/dis_P(mx, my))*exp(-2*1j*pi*(dis_P(mx, my)-dis_0(k)) / lamda);
        end
    end
    
    H2(:, k) = reshape(Arr_U',1, []);

end


end

