function [ur,ut,up,dt,npts] =syn_modes_funcM_comp(modesS,...
    modesT,rs,thetas,phis,rr,thetar,phir,M,f,dis,vp,vs,Q,dt,npts)
dd=1;
xr=rr*sin(thetar)*cos(phir);
yr=rr*sin(thetar)*sin(phir);
zr=rr*cos(thetar);
f1=M/dd*f;
[ff]=force_convert(f1,thetas,phis);
[phir1,thetar1,rr1]=cart2sph(xr+dd/2*dis(1),yr+dd/2*dis(2),zr+dd/2*dis(3));
thetar1=pi/2-thetar1;
[ur1,ut1,up1]=syn_modes_func3(modesS,modesT,rs,thetas,...
    phis,rr1,thetar1,phir1,ff,vp,vs,Q,dt,npts);

[phir2,thetar2,rr2]=cart2sph(xr-dd/2*dis(1),yr-dd/2*dis(2),zr-dd/2*dis(3));
thetar2=pi/2-thetar2;
[ur2,ut2,up2]=syn_modes_func3(modesS,modesT,rs,thetas,...
    phis,rr2,thetar2,phir2,-ff,vp,vs,Q,dt,npts);
ur=ur1+ur2;
ut=ut1+ut2;
up=up1+up2;
end