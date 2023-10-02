%% ATEHE for step size and time
%The effect of step size on time and NMSE
clc;clear;
%close all;

%%  %%%%%%%%%%%%%%%%%%%Parameters setting%%%%%%%%%%%%%%%%%%%%%%
K=8;                        %number of user
L=64;                       %equal channel length

m = 128;                     %number of received antennas
KL = K*L;                    
Np = 256;                    %length of pilots

KB = 64;                     %number of non-zero rows
p = 45/64;                   %in-row sparsity ratio
SNR = [15,30];                %signal noise ratio

lsnr = length(SNR);

itermax=10000;               
prin = 1000;                
num_A = 6;                 

s = 0:8:6*K;                  %testing step size
s(1)=1;
len = length(s);
all_time = zeros(1,lsnr);

NMSE = zeros(lsnr,len);
Time = zeros(lsnr,len);

%% %%%%%%%%%%%%%%%%%%%Starting Simulation%%%%%%%%%%%%%%%%%%%%%%
tic;
for snr_index=1:lsnr
    nmse = zeros(1,len);

    time = zeros(1,len);
    temp = zeros(1,len);
    SNR_tt = SNR(snr_index);
for tt=1:itermax
%！！！！！！！！！！！！！！！！！！！！ Generating data！！！！！！！！！！！！！！！！！！！！！！！！
        time_temp =zeros(1,len);
        nmse_temp = zeros(1,len);
        [Y,S,X,W,Index_KB] = Gen_Data(K,L,m,Np,p,KB,SNR_tt);

     for ins = 1:len
        tic;   
        [H_r,ZSN] = ATEHE( Y,S,s(ins),SNR_tt,W,Index_KB,X);            %6A-TEHE
        nmse_temp(ins) = norm(H_r-X,'fro')^2/norm(X,'fro')^2;
        
        time_temp(ins) = toc;
     end
     nmse = nmse+nmse_temp;
     time = time + time_temp;

     txt=['tt=',num2str(tt), '  nmmse:',num2str(nmse_temp)]; 

     if rem(tt,prin) == 0                                                    
      disp(txt);
     end

end
all_time(snr_index) = sum(time);
Time(snr_index,:) = time /itermax;
NMSE(snr_index,:) = nmse/itermax;
end

figure('Name','Step-Time-NMSE3D')
plot3(s,Time(1,:),NMSE(1,:),'.r-');
hold on;
plot3(s,Time(2,:),NMSE(2,:),'+b-');
set(gca,'zscale','log');
grid on;
xlabel('Stepsize');
ylabel('Time');
zlabel('NMSE');

%==========2D Time vs stepsize====================
figure('Name','Step-Time2D')
plot(s,Time(1,:),'.r-');
hold on;grid on;
plot(s,Time(2,:),'+b-');
xlabel('Stepsize');
ylabel('Time');

%save l1stn1.mat;
