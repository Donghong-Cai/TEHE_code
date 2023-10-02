%% Threshold-Enhanced Hierarchical Spatial Non-Stationary Channel Estimation
%Y=XH+W£¬Known X£¬Y  estimating H£¬Y---Np*m, S---Np*KL, X---KL*m
clear ; clc;
close all;

%%  %%%%%%%%%%%%%%%%%%%Parameters setting%%%%%%%%%%%%%%%%%%%%%%
K=8;                        %number of user
L=64;                       % equal channel length

m = 128;                     % number of received antennas
KL = K*L;                    
Np = 256;                    %length of pilots

KB = 64;                     %number of non-zero rows
p = 45/64;                   %in-row sparsity ratio
SNR = 0:5:30;                %signal noise ratio

lsnr = length(SNR);

itermax=10000;               
prin = 200;                 %cycle for output
num_A = 5;                  %number of algorithms
NMSE = zeros(lsnr,num_A);           
supp_rate = zeros(lsnr,num_A);
supp_all = zeros(lsnr,num_A);
supp_11 = zeros(lsnr,num_A);
supp_00 = zeros(lsnr,num_A);

ZSNF = zeros(lsnr,num_A);        

%% %%%%%%%%%%%%%%%%%%%Starting Simulation%%%%%%%%%%%%%%%%%%%%%%
tic;
for snr_index=1:lsnr
    nmse=zeros(1,num_A);
    temp_ZN=zeros(1,num_A);
    snr=10^(-SNR(snr_index)/10);
    SNR_tt = SNR(snr_index);

    for tt=1:itermax
%------------------------------------------- Generating data-------------------------------------------
        [Y,S,X,W,Index_KB] = Gen_Data(K,L,m,Np,p,KB,SNR_tt);
        ZSN = zeros(1,num_A);
 %%   -------------------------------------------Channel Estimation-------------------------------------------
        H_r = zeros(KL,m,num_A);

        [H_r(:,:,1),H_r(:,:,2),ZSN(2)] = TEHE( Y,S,1,SNR_tt,W,Index_KB,X);     %1:Block-SAMP£¬2:TEHE
      
        for ii=1:m
            H_r(:,ii,3) = Parallel_SAMP( Y(:,ii),S,1,SNR_tt);                  %3Parallel-SAMP   -large complexity, individual testing better
        end

        for ii=1:m
           H_r(:,ii,4) = oracle_Ls( Y(:,ii),S,X(:,ii));                       %4oracle_LS
        end

        [H_r(:,:,5),ZSN(5)] = ATEHE( Y,S,K,SNR_tt,W,Index_KB,X);            %5A-TEHE


   %% -------------------------------------------Data Analysis-------------------------------------------

        H1=abs(X);
        %-------------------------------------------debug-------------------------------------------
        temp_ZN=temp_ZN+ZSN;
        %-------------------------------------------success rate for recovery-------------------------------------------
        supp_rate_tem = zeros(lsnr,num_A);
        supp_all_tem = zeros(lsnr,num_A);
        supp_11_tem = zeros(lsnr,num_A);
        supp_00_tem = zeros(lsnr,num_A);

        for ii=1:num_A
            Ha_r=abs(H_r(:,:,ii));
            temp= sum(sum(~xor(H1(Index_KB,:),abs(Ha_r(Index_KB,:)))));
            supp_rate_tem(snr_index,ii)=supp_rate_tem(snr_index,ii)+temp;
            supp_all_tem(snr_index,ii)=supp_all_tem(snr_index,ii)+sum(sum(~xor(H1,Ha_r)));    
            supp_11_tem(snr_index,ii)=supp_11_tem(snr_index,ii)+length(Ha_r(Ha_r(H1>0)>0)); %1-->0
            supp_00_tem(snr_index,ii)=supp_00_tem(snr_index,ii)+length(Ha_r(Ha_r(H1==0)==0)); %0-->1

        end
        supp_rate = supp_rate+supp_rate_tem;
        supp_all = supp_all+supp_all_tem;
        supp_11 = supp_11+supp_11_tem;
        supp_00 = supp_00+supp_00_tem;
        %------------------------------NMSE--------------------------------
        nmse_1=zeros(1,num_A);
        for jj=1:num_A
            nmse_1(1,jj) = sum(sum(abs((H_r(:,:,jj)-X)).^2))/sum(sum(abs(X).^2));  
        end
        nmse = nmse+nmse_1;                                                       
        txt=['t=',num2str(tt),' nmse:',num2str(nmse_1),' total:',num2str(nmse)];
        if rem(tt,prin)==0                                                          %display
            disp(txt);
        end

    end
    toc
    ZSNF(snr_index,:)=temp_ZN;
    NMSE(snr_index,:)=nmse./itermax;                       
end

%disp(NMSE);                                  
supp_rate = supp_rate/(itermax*KB*m);
supp_all = supp_all/(itermax*KL*m);
supp_11 = supp_11/(itermax*KB*m*p);
supp_00 = supp_00/(itermax*(KL*m-KB*m*p));

%% %%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%
%===========SNRvsNMSE all algorthim====================
figure('Name','NMSE-all');
semilogy(SNR,NMSE(:,3),'.k-', SNR,NMSE(:,2),'xm-', SNR,NMSE(:,4),'+r--', SNR,NMSE(:,5),'db-')

grid on;

xlabel('SNR(dB)');
ylabel('NMSE');
legend('Parallel-SAMP','Proposed TEHE','Proposed A-TEHE','Orcale-LS');

%===========Success rate for recovery ====================
figure('Name','SUPPORT');
plot(SNR,supp_rate(:,3),'.k-',SNR,supp_rate(:,2),'xm-', SNR,supp_rate(:,4),'+r--', SNR,supp_rate(:,5),'db-');
grid on;
xlabel('Success rate for recovery');
ylabel('SUPPORT');

legend('Parallel-SAMP','Proposed TEHE','Proposed A-TEHE','Orcale-LS');

%===========SNRvsNMSE second layer====================
figure('Name','2ndlayer-NMSE'); 
semilogy( SNR,NMSE(:,1),'*g-', SNR,NMSE(:,2),'xm-', SNR,NMSE(:,5),'db-', SNR,NMSE(:,4),'+r--');
grid on;

xlabel('SNR(dB)');
ylabel('NMSE');

legend('Block-SAMP','Proposed TEHE','Proposed A-TEHE','Orcale-LS');

%===========Success rate for recovery for second layer --non-zero elements ====================
figure('Name','SUPPORT_11');
plot(SNR,supp_11(:,1),'*g-',SNR,supp_11(:,2),'xm-', SNR,supp_11(:,5),'db-');
grid on;
xlabel('SNR(dB)')
ylabel('Success rate for recovery')
legend('Block-SAMP','Proposed TEHE','Proposed A-TEHE');

%===========Success rate for recovery for second layer --zero elements ====================
figure('Name','SUPPORT_00');
plot(SNR,supp_00(:,1),'*g-',SNR,supp_00(:,2),'xm-', SNR,supp_00(:,5),'db-');
grid on;
xlabel('SNR(dB)')
ylabel('Success rate for recovery')
legend('Block-SAMP','Proposed TEHE','Proposed A-TEHE');

%% Saving
%   save SAMP.mat
