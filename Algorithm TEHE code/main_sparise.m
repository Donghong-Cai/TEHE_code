%% second layer analysis
% Algorithm accuracy versus sparsity within non-zero rows
clear ; clc;
close all;

%%  %%%%%%%%%%%%%%%%%%%Parameters setting%%%%%%%%%%%%%%%%%%%%%%
K=8;                        %number of user
L=64;                       %equal channel length

m = 128;                     %number of received antennas
KL = K*L;                    
Np = 256;                    %length of pilots

KB = 64;                     %number of non-zero rows
p = (16:8:56)/64;                   %in-row sparsity ratio
SNR = 30;                %signal noise ratio

l_p = length(p);

itermax=10000;               
prin = 1000;                 %cycle for output
num_A = 5;                  %number of algorithms
NMSE = zeros(l_p,num_A);            

ZSNF = zeros(l_p,num_A);        
%% %%%%%%%%%%%%%%%%%%%Starting Simulation%%%%%%%%%%%%%%%%%%%%%%
tic;
for p_index=1:l_p
    nmse=zeros(1,num_A);
    temp_ZN=zeros(1,num_A);
    p_t = p(p_index);
    SNR_tt = SNR;

    for tt=1:itermax
%-------------------------Generating Data---------------------------
        [Y,S,X,W,Index_KB] = Gen_Data(K,L,m,Np,p_t,KB,SNR_tt);
        ZSN = zeros(1,num_A);
  %%   -------------------------Channel Estimation-------------------------------------
        H_r = zeros(KL,m,num_A);
 
         [H_r(:,:,1),H_r(:,:,2),ZSN(2),lambda] = TEHE( Y,S,1,SNR_tt,W,Index_KB,X);     %1:Block-SAMPï¼?2:TEHE
      
        for ii=1:m
            H_r(:,ii,3) = Parallel_SAMP( Y(:,ii),S,1,SNR_tt);                  %3Parallel-SAMP   -large complexity, individual testing better 
        end

        for ii=1:m
           H_r(:,ii,4) = oracle_Ls( Y(:,ii),S,X(:,ii));                       %4oracle_LS
        end
% 
         [H_r(:,:,5),ZSN(5)] = ATEHE( Y,S,K,SNR_tt,W,Index_KB,X);            %5A-TEHE


   %% --------------------------------data analysis--------------------------------
        %-------------------NMSE----------------
        nmse_1=zeros(1,num_A);
        for jj=1:num_A
            nmse_1(1,jj) = sum(sum(abs((H_r(:,:,jj)-X)).^2))/sum(sum(abs(X).^2));  
        end
        nmse = nmse+nmse_1;                                                        
        %txt=['t=',num2str(tt),' nmse:',num2str(nmse_1),' total:',num2str(nmse)];
        txt=['t=',num2str(tt),' nmmse:',num2str(nmse_1)];
        if rem(tt,prin)==0                                                          
            disp(txt);
        end

    end
    toc
    ZSNF(p_index,:)=temp_ZN;
    NMSE(p_index,:)=nmse./itermax;                        
end

%% %%%%%%%%%%%%%%%%%%%%%%%plot NMSE vs. in-row%%%%%%%%%%%%%%%%%%%%%%%%%%
p = p*m;
figure('Name','NMSE vs sparisity of non-zero rows');
semilogy(p,NMSE(:,1),'*b-',p,NMSE(:,2),'xm-',p,NMSE(:,5),'dc-',p,NMSE(:,4),'+r--')
hold on; grid on;
xlabel('sparisity of non-zero rows')
ylabel('NMSE')
legend('Block-SAMP','Proposed TEHE','Proposed A-TEHE','Orcale-LS');

%% saving
 % save sprase.mat
