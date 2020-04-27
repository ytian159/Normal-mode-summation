function [Ur,Vr]=RadSphMod2(l,ka,kb,r,x1,x2)
Ur=ka.^(-1).*kb.^(-1).*r.^(-1).*(kb.*l.*SphericalBesselJ(l,ka.*r)+ ...
  ka.*(l.*(1+l).*x1.*SphericalBesselJ(l,kb.*r)+(-1).*kb.*r.* ...
  SphericalBesselJ(1+l,ka.*r)));
Vr=ka.^(-1).*kb.^(-1).*(l.*(1+l)).^(1/2).*r.^(-1).*(kb.* ...
  SphericalBesselJ(l,ka.*r)+ka.*x2.*((1+l).*SphericalBesselJ(l,kb.*r) ...
  +(-1).*kb.*r.*SphericalBesselJ(1+l,kb.*r)));
% nmmod=(Ur.^2+Vr.^2).*r.^2*rho*dr;
% nmfac=sqrt(1/sum(nmmod));
% UU=Ur*nmfac;
% VV=Vr*nmfac;


end
