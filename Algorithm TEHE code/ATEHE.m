function [ Theta1, ZSN ] = ATEHE( Y,X,S,SNR,w,supp,H)
%A-TEHE
%ZSN---debug

    [~,m] = size(Y);
    [t,n] = size(X);        %sensing matrix
    Theta1 = zeros(n,m);    %saving theta(columns vector)
    Pos_theta = [];         %index of slected columns
    R_n = Y;                %initializing residual R_0=y
    L = S;                  %Size of the finalist in the first stage
    Stage = 1;              %initializing Stage
    IterMax = t;

    p_w=t*m/10^(SNR/10);         %power of noise
    epsilon = p_w;
    for ii=1:floor(IterMax/S)
        Product = X'*R_n;                               %The inner product of the columns of the sensing matrix A with the residuals
        [~,pos] = sort(sum(abs(Product),2),'descend'); 
        Sk = pos(1:L);                                  %The most L large
        Ck = union(Pos_theta,Sk);
        
        if length(Ck)<=t
            At = X(:,Ck);                               
        else
            Theta_ls=0;
            break;
        end

        Theta_ls = (At'*At)^(-1)*At'*Y;   %LS solving
        [~,pos]=sort(sqrt(sum(Theta_ls.^2,2)),'descend');
        
        F = Ck(pos(1:L));
        Theta_ls = (X(:,F)'*X(:,F))^(-1)*X(:,F)'*Y;
        R_new = Y - X(:,F)*Theta_ls;            %updating residual

         X_s=X(:,F);

        if norm(R_new,"fro")^2 < epsilon    %halting condition true 
            Pos_theta = F;
            break;
        elseif norm(R_new,"fro")>=norm(R_n,"fro")
            Stage = Stage + 1;
            L = Stage*S;
            if ii == IterMax
                Pos_theta = F;
            end
            ii = ii - 1;
        else
            Pos_theta = F;
            R_n = R_new;
        end
    end
      ZS=setdiff(supp,Pos_theta);                %omitted support set 
      ZSN=length(ZS);  
      s = length(F);
%-------------enhance layer------------------------
      lambda=norm((X_s'*X_s)^(-1)*X_s')/sqrt(10^(SNR/10)*s/t)*sqrt(SNR/30);
      Theta1(F,:)=Theta_ls;
      Theta1(abs(Theta1)<lambda)=0;

end