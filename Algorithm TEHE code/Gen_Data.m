function [Y,S,X,W,Index_KB] = Gen_Data(K,L,M,Np,p,KB,SNR)
%Y--resive signal, S--sencing matrix, X--channel matrix, W--noise matrix
%K--number of users, L--length of OFDM, M--resive anntant, Np--length of pilot
%p--sprasity of non-zero rows, KB--number of non-zero rows, snr--signal noise radio
%This channel is randomly generated based on the characteristic of inhomogeneous sparsity, 
%which is more mathematically general. 
%Other analog channels with this characteristic can also be used

%！！！！！！！！！！！！！！！！！！！！ Generating Data！！！！！！！！！！！！！！！！！！！！！！！！
   % ！！！！！！Channel and noise！！！！
        KL = K*L;
        X=zeros(KL,M);
        hpre = zeros(1,KB*M);
        Index_KB = randperm(KL,KB);                             %select KB non-zero rows
        h_index = randperm(KB*M,p*KB*M);                        
        noz =(randn(1,p*KB*M)+randn(1,p*KB*M)*1i);              %random non-zero element
        hpre(1,h_index) = noz;                                  %insert
        hpre=reshape(hpre,KB,M);                                
        X(Index_KB,:)=hpre;                                     
        
        snr = 10^(-SNR/10);
        W = sqrt(snr/2)*(randn(Np,M)+randn(Np,M)*1i);           %awang

   %！！！！！！！！！！Pilots！！！！！！！！！！！！！！！！
        
       %--------------Gauss-----------------
           S = randn(Np,KL);
       %-------------Fourier------------ 
%          S = zeros(Np,KL);
%         for ii = 1:K
%             Theta = rand(Np,1)*pi*2;
%             pk = exp(1i*Theta);
%             Pk =diag(pk);
%             FFt = fft(eye(Np));
%             S(:,(ii-1)*L+1:ii*L) = Pk*FFt(:,1:L);
%         end

%           Theta = rand(KL,1)*pi*2;
%           pk = exp(1i*Theta);
%           Pk = diag(pk);
%           FFt = fft(Pk);
%           S = FFt(randi(KL,Np,1),:);

       %-----------ZC pilot-----------------
%         r = randi(Np-1,1,K);
%         for ii = 1:K
%             pk = zadoffChuSeq(r(ii),Np);
%             Pk =diag(pk);
%             FFt = fft(eye(Np));
%             S(:,(ii-1)*L+1:ii*L) = Pk*FFt(:,1:L);
%         end
         for ii=1:KL                        %normalization
            S(:,ii) = S(:,ii)/norm(S(:,ii));
         end
%           tes = abs(S(:,1)'*S(:,500));
        %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
        % rn=norm(S,"fro")^2;
        % rw=norm(W,"fro")^2;
        %X = X*diag(1./sqrt(sum(X.*X)));
        Y = S*X+W; 
end