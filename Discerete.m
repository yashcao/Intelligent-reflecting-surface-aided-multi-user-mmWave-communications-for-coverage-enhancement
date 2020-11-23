function [ Theta_dis] = Discerete( bit, Theta_opt )
%DISCERETE 此处显示有关此函数的摘要

F_set = 2*pi*(0:(2^bit-1))/(2^bit);

Theta_ang = angle(Theta_opt)+pi;

% len = length(Theta_ang);
% Theta_dis = zeros(len, 1);
% 
% for i=1:len
%     Theta_d = abs(Theta_ang(i)-F_set);
%     [~, I] = min(Theta_d);
%     Theta_dis(i) = exp(1i*F_set(I));
% end

size_Theta = size(Theta_ang);
Theta_dis = zeros(size_Theta);

for i=1:size_Theta(1)
    for j=1:size_Theta(2)
        Theta_d = abs(Theta_ang(i,j)-F_set);
        [~, I] = min(Theta_d);
        Theta_dis(i,j) = exp(1i*F_set(I));
    end
end


end

