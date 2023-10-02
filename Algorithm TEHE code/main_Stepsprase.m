%% ATEHE for step size and sparsity
%The relationship among step size, row sparsity and NMSE
clc;clear;
close all;

%%  %%%%%%%%%%%%%%%%%%%Parameters setting%%%%%%%%%%%%%%%%%%%%%%
K=8;                        %number of user
L=64;                       %equal channel length

m = 128;                     %number of received antennas
KL = K*L;                    
Np = 256;                    %length of pilots

KB = 24:8:80;                     %testing number of non-zero rows
lenk = length(KB);
p = 45/64;                   %in-row sparsity ratio
SNR = [15,30];                %signal noise ratio

lsnr = length(SNR);

itermax=1;               
prin = 1000;                 


s = 0:8:10*K;                  %testing step size
s(1)=1;
len = length(s);
nmse = zeros(lenk,len);

%% %%%%%%%%%%%%%%%%%%%Starting Simulation%%%%%%%%%%%%%%%%%%%%%%
for snr_index=1:lsnr
    SNR_tt = SNR(snr_index);
    tic;
    for ikb = 1:lenk
        nmse_temp = zeros(1,len);
        KB_tt =  KB(ikb);

        for tt=1:itermax
            %！！！！！！！！！！！！！！！！！！！！  Generating data！！！！！！！！！！！！！！！！！！！！！！！！
            [Y,S,X,W,Index_KB] = Gen_Data(K,L,m,Np,p,KB_tt,SNR_tt);
            temp = zeros(1,len);
            for ins = 1:len
                [H_r,ZSN] = ATEHE( Y,S,s(ins),SNR_tt,W,Index_KB,X);            %A-TEHE
                temp(ins) = norm(H_r-X,'fro')^2/norm(X,'fro')^2;
            end
            nmse_temp = nmse_temp + temp;
            %nmseT = nmseT + nmse_temp;
            %txt=['tt=',num2str(tt),' KB=', num2str(KB(ikb)), ' nmmsetem:',num2str(nmse_temp),' nmse=',num2str(temp)];
            if rem(tt,prin)==0
                txt=['tt=',num2str(tt),' ||X||=',num2str(KB_tt),' nmse=',num2str(temp)];
                disp(txt); 
            end
        end

        nmse(ikb,:)  = nmse_temp/itermax;
        txt=[' ||X||=', num2str(KB_tt), ' nmmse:',num2str(nmse(ikb,:))] ;
        disp(txt); 
        toc
    end
    %=============saving====================
    if SNR_tt == 15
       save l1sxn1w15.mat;
    else
       save l1sxn1w30.mat;
    end

     %=============plot step size, row sparsity VS NMSE===================
    [A,B] = meshgrid(KB,s);
    figure('Name','Step-||X||_{2,0}-NMSE')
    surf(A,B,nmse');
    set(gca,'zscale','log');
    grid on;
    xlabel('KB');
    ylabel('Step');
    zlabel('NMSE');
end

