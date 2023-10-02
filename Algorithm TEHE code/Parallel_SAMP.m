function [ theta ] = Parallel_SAMP( y,A,S,rio )
%Parallel-SAMP--Estimation using SAMP for each antenna individually

    [M,N] = size(A);            
    [t,~] = size(y);
    theta = zeros(N,1);         
    Pos_theta = [];             
    r_n = y;                    
    L = S;                      
    Stage = 1;                  
    IterMax = M;

    p_w = t/10^(rio/10);
    epsilon = p_w;								
    for ii=1:IterMax                                    
        product = A'*r_n;                               
        [~,pos]=sort(abs(product),'descend');           
        Sk = pos(1:L);                                  
        Ck = union(Pos_theta,Sk);                      

        if length(Ck)<=M
            At = A(:,Ck);                               
        else
            theta_ls=0;
            break;
        end

        theta_ls = (At'*At)^(-1)*At'*y;                 
        [~,pos]=sort(abs(theta_ls),'descend');         
        F = Ck(pos(1:L));                               
        
        X_s = A(:,F);
        theta_ls = (X_s'*X_s)^(-1)*X_s'*y;
        r_new = y - A(:,F)*theta_ls;                    
        
       
        if norm(r_new)^2 < epsilon                      
            Pos_theta = F;
            break;
        elseif norm(r_new) >= norm(r_n)                   
            Stage = Stage + 1;
            L = Stage*S;
            if ii == IterMax                
                Pos_theta = F;              
            end
            ii = ii - 1;
        else
            Pos_theta = F;
            r_n = r_new;
        end
    end
    theta(Pos_theta)=theta_ls;
end