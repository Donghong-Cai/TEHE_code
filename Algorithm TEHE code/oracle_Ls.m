function [theta] = oracle_Ls(y, S, x)
%Optimal Ls estimation with known correct support set
[KL,~] = size(x);
theta = zeros(KL,1);
supp = (x~=0);
S_e = S(:,supp);  
theta(supp) =  (S_e'*S_e)^(-1)*S_e'*y;

end