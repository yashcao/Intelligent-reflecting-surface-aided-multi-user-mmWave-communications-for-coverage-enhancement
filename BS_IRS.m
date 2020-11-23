function [G] = BS_IRS(Mx,My,Nx,Ny,L,d)
% BS-IRS channel matrix
N = Nx*Ny;
M  = Mx*My;
%G = zeros(M, N);

lambda_gain = 9.60; % 9.82 dBi
a = 61.4;
b = 2;
sigma = 3.8;  % 5.8 dB
kappa = a+10*b*log10(d)+normrnd(0, sigma, 1);
% channel attenuation
var = [1, 0.01*ones(1, L)].* (sqrt(1/2*10^(-0.1*kappa)))*(randn(1)+1j*randn(1));



% ULA array steering vector
theta = rand(1, L+1)-0.5;
Arr_BS = exp(-2*1i*pi*theta'*((0:N-1)-(N-1)/2))';
%Arr_BS = exp(-2*1i*pi*theta'*(0:N-1))';


%beta_BL = sqrt(var/2).*(randn(1, L+1) +1i*randn(1, L+1));
beta_BL = lambda_gain*var;

% random distribution of LIS angles parameters
angL_e = rand(1, L+1)-0.5;                          % elevation sin_e
%angL_a = (rand(1, L+1)-0.5).*cos(asin(angL_e));    % azimuth cos_a*cos_e
angL_a = rand(1, L+1)-0.5;

% array steering vector
arrL_a = (1/sqrt(Mx))*exp(-2*1i*pi*angL_a'*((0:Mx-1)-(Mx-1)/2))';
arrL_e = (1/sqrt(My))*exp(-2*1i*pi*angL_e'*((0:My-1)-(My-1)/2))';


for l=1:L+1
    Arr_LIS(:,l) = kron(arrL_a(:,l), arrL_e(:,l));
end

Beta = diag(beta_BL);
G = sqrt(M*N)*Arr_LIS*Beta*Arr_BS';

end

