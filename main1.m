% UPA_reflecting_surface
clear;close all;clc;
%% BS parameters settings
Nx = 8;    % TX antenna nums
Ny = 4;    % TX antenna nums
K = 2;     % user nums
L = 2;     % NLOS path nums
%SNR = 30; % SNR in dB
N = Nx*Ny;

d_BI = 50;
cord = [(2*rand(K,1)-1)+4, 2*rand(K,1)-1]*10;
for ii=1:K
    c_IU(ii,:) = cord(ii,:)-[40, 30];
end

d_IU = sqrt(sum(c_IU.^2, 2));


%% surface parameters
Mx = 5;                      % x-axis surface element nums
My = [1,2,3,4,5,6]*2;   % y-axis surface element nums
M  = Mx*My;


%% optimization parameters
P_max = 1; % W (30dBm)
sigma = 3.16e-12; % dBm-->W (-85dBm)

SR = zeros(length(My),1);
SR1 = zeros(length(My),1);
SR2 = zeros(length(My),1);
SR_rand = zeros(length(My),1);

%% plot
for mi=1:length(My)
for test = 1:1000

%% BS-IRS channel matrix
G = BS_IRS(Mx,My(mi),Nx,Ny,L,d_BI);

%% IRS-user channel matrix
H = FF_User(Mx,My(mi),K,d_IU);

%% AO algorithm
[ Sum_Rate, P, Theta ] = Opt_func( M(mi),N,K,P_max,sigma,H,G );

DisSumRate1 = DiscreteSum( 1,Theta,K,P,sigma,H,G );
DisSumRate2 = DiscreteSum( 2,Theta,K,P,sigma,H,G );

Sum_rand = RandSum( 0, M(mi),N,K,P_max,sigma,H,G );


SR(mi) = SR(mi)+Sum_Rate;
SR1(mi) = SR1(mi)+DisSumRate1;
SR2(mi) = SR2(mi)+DisSumRate2;
SR_rand(mi) = SR_rand(mi)+Sum_rand;
end

SR(mi) = SR(mi)/test;
SR1(mi) = SR1(mi)/test;
SR2(mi) = SR2(mi)/test;
SR_rand(mi) = SR_rand(mi)/test;
end


figure(1);
plot(M, SR, 'b-.o', 'LineWidth',1,'MarkerSize',8);
hold on;
plot(M, SR1, '--s', 'color', [0.13 0.55 0.13], 'LineWidth',1,'MarkerSize',10);
plot(M, SR2, 'r--^', 'LineWidth',1,'MarkerSize',8);
plot(M, SR_rand, 'm-*', 'LineWidth',1,'MarkerSize',10);
set(gca,'xtick',M);
xlim([M(1), M(end)]);
xlabel('M'),ylabel('Sum-rate'); %xlim([1, iter]);
grid on;
handle=legend('Proposed, $\mathcal{F}_{\rm c}$','Proposed, $\mathcal{F}_{\rm d}$, $b=1$',...
    'Proposed, $\mathcal{F}_{\rm d}$, $b=2$', 'ZF+Random PBF');
set(handle, 'Interpreter', 'LaTex');
set(handle, 'FontSize', 13, 'Location', 'best');
set(gca, 'FontSize', 12);



