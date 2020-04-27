function [Sr,St,Sp]=Slm_func(l,m,theta,phi)
%% calculate the three components of Spherical Harmonics Slm

Sr=0;
St=(l.*(1+l)).^(-1/2).*(m.*cot(theta).*ylm(l,m,theta,phi)+exp(1).^(( ...
  sqrt(-1)*(-1)).*phi).*(1+l+(-1).*m).^(1/2).*(2+l+m).^(1/2).*ylm(l,1+m, ...
  theta,phi));
Sp=sqrt(-1).*(l.*(1+l)).^(-1/2).*m.*csc(theta).*ylm(l,m,theta,phi);
if l==0
    Sr=0;
    St=0;
    Sp=0;
end
% if theta==0
%     Sr=0;
%     St=asso_P_theta0_coef(l,m)*theta^(m-1)*m;
% end
end