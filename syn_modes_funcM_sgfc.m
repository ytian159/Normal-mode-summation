function [ur,ut,up,dt,npts] =syn_modes_funcM_sgfc(modesS,...
    modesT,rs,thetas,phis,rr,thetar,phir,M,vp,vs,Q,dt,npts)
dd=1;
xr=rr*sin(thetar)*cos(phir);
yr=rr*sin(thetar)*sin(phir);
zr=rr*cos(thetar);
ur=0;
ut=0;
up=0;

f1=M(1,1)/dd*[1;0;0];
[phir1,thetar1,rr1]=cart2sph(xr+dd/2,yr,zr);
thetar1=pi/2-thetar1;
[usynr,usynt,usynp]=syn_modes_func3(modesS,modesT,rs,thetas,...
    phis,rr1,thetar1,phir1,f1,vp,vs,Q,dt,npts);
ur=ur+usynr;
ut=ut+usynt;
up=up+usynp;
f2=M(1,1)/dd*[-1;0;0];
[phir2,thetar2,rr2]=cart2sph(xr-dd/2,yr,zr);
thetar2=pi/2-thetar2;
[usynr,usynt,usynp]=syn_modes_func3(modesS,modesT,rs,thetas,...
    phis,rr2,thetar2,phir2,f2,vp,vs,Q,dt,npts);
ur=ur+usynr;
ut=ut+usynt;
up=up+usynp;
end