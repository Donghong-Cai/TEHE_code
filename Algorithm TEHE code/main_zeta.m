%% kappa relationship
%Y=XH+W£¬Known X£¬Y  estimating H£¬Y---Np*m, S---Np*KL, X---KL*m
clear ; clc;
close all;

%%  %%%%%%%%%%%%%%%%%%%Parameters setting%%%%%%%%%%%%%%%%%%%%%%
K=8;                        %number of user
L=64;                       %equal channel length

m = 128;                     %number of received antennas
KL = K*L;                    
Np = 256;                    %length of pilots

KB = 64;                     %number of non-zero rows
p = 45/64;                   %in-row sparsity ratio
SNR = 0:5:30;                %signal noise ratio

lsnr = length(SNR);

itermax=10000;               

c = 0:0.05:1;						%tested kappa
num_L = length(c);                  %number of tested kappa
NMSE = zeros(lsnr,num_L);            
supp_rate = zeros(lsnr,num_L);
supp_all = zeros(lsnr,num_L);
ZSNF = zeros(lsnr,num_L);       


%% %%%%%%%%%%%%%%%%%%%Starting Simulation%%%%%%%%%%%%%%%%%%%%%%
tic;
for snr_index=1:lsnr
    nmse=zeros(1,num_L);
    temp_ZN=zeros(1,num_L);
    snr=10^(-SNR(snr_index)/10);
    SNR_tt = SNR(snr_index);

    for tt=1:itermax
        %--------------------------------------- Generating Data-------------------------------------------
        [Y,S,X,W,Index_KB] = Gen_Data(K,L,m,Np,p,KB,SNR_tt);
        ZSN = zeros(1,num_L);
        %----------------------------------------Channel Estimation-----------------------------------------?
        [X_r1,X_r2,ZSN(2),lambda] = TEHE( Y,S,1,SNR_tt,W,Index_KB,X);     

        H1=abs(X);

        % ---------------------------kappa(lambda)-------------------------------------------
        XH = zeros(KL,m,num_L);
        for ii = 1:num_L
            BM = X_r1;
            zeta = lambda*c(ii);
            BM(abs(BM)<zeta)=0;
            XH(:,:,ii) = BM;
            nmse(ii)=nmse(ii)+norm(BM-X,"fro")^2/norm(X,"fro")^2;
        end
    end
    toc
    NMSE(snr_index,:)=nmse./itermax;
end
save kappa.mat;

[A,B] = meshgrid(SNR,c);
figure('Name','NMSE-C_SNR')
surf(A,B,NMSE');
set(gca,'zscale','log');
grid on;
xlabel('SNR');
ylabel('Parameter $\kappa$');
zlabel('NMSE');