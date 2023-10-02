function [ Theta,ThetaH, ZSN,lambda] = TEHE( Y,X,S,SNR,w,supp,H)
%Theta---Block-SAMP£¬ThetaH---TEHE£¬
%ZSN----debug,lambda--threshold

[~,m] = size(Y);
[t,n] = size(X);            %sensing matrix
Theta = zeros(n,m);         %saving theta(columns vector)
Pos_theta = [];             %index of slected columns
R_n = Y;                    %initializing residual R_0=y
L = S;                      %Size of the finalist in the first stage
Stage = 1;                  %initializing Stage
IterMax = t;                

p_w = t*m/10^(SNR/10);        %power of noise
epsilon = p_w;
for ii = 1:IterMax            
    Product = X'*R_n;                                           %The inner product of the columns of the sensing matrix A with the residuals
    [~,pos] = sort(sum(abs(Product),2),'descend');              
    Sk = pos(1:L);                                              %%The most L large
    Ck = union(Pos_theta,Sk);                                  
    
    if length(Ck)<=t
        At = X(:,Ck);                                           
    else
        Theta_ls = 0;
        break;
    end

    Theta_ls = (At'*At)^(-1)*At'*Y;                     %LS solving
    [~,pos] = sort(sqrt(sum(Theta_ls.^2,2)),'descend');   
    F = Ck(pos(1:L));
    X_s = X(:,F);
    Theta_ls = (X_s'*X_s)^(-1)*X_s'*Y;
    R_new = Y - X(:,F)*Theta_ls;                        %%updating residual
    
    if norm(R_new,"fro")^2 < epsilon                              %halting condition true
        Pos_theta = F;
        break;
    elseif norm(R_new,"fro") >= norm(R_n,"fro")                   
        Stage = Stage + 1;                                        %Update the stage index
        L = Stage*S;                                              %Update the size of finalist
        if ii == IterMax                                          
            Pos_theta = F;                                        
        end
        ii = ii - 1;
    else
        Pos_theta = F;                  %Update the finalist Fk
        R_n = R_new;                    %Update the residue
    end
end
ZSN = length(setdiff(supp,F));  

Theta(Pos_theta,:) = Theta_ls;

%-------------enhance layer------------------------
ThetaH = Theta;
s = length(F);
invx=(X_s'*X_s)^(-1)*X_s';

 lambda=norm(invx)/sqrt(10^(SNR/10)*s/t); 
 ThetaH(abs(ThetaH)<lambda*sqrt(SNR/30))=0;
end