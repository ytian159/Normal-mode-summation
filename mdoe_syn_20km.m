%%
clearvars;
set(0,'defaultaxesfontsize',16); % control the default fontsize of the plots
set(0,'defaulttextfontsize',16);
%load('modes_data_20km_3hz.mat');
load('modes_data_20km_1p5hz.mat');
vp = 6e3;
vs = 3e3;

rho = 3000;
mu = rho* vs^2;
lamda = rho*vp^2 - 2*mu;
R= 20e3;
Q=50;
T=50;
%%
f0=0.5;
dt=0.1;
fs=1/dt;
nt = round(T/dt);
%nt = nt + mod(nt,2);
time = [0:nt-1]*dt;
df = 1/(nt*dt);
dw = 2*pi*df;

fmax = 3*f0; % hz

ts =1.5/f0; 
 [w,t] = ricker(f0,ts,dt,nt); 
 wf = ifft(w)*dt*nt; 
 %%
M=[1,0,0;
   0,1,0;
   0,0,1];
fsrc=[1,1,1];
%R=6371e3;
Rc=R/2;
h=R*0.5;
rs = R-h; % m radius
thetas = 0.5*pi;
phis = 1*pi;
xs=rs*sin(thetas)*cos(phis);
ys=rs*sin(thetas)*sin(phis);
zs=rs*cos(thetas);

ntr=36;
phiar=(1:ntr)/ntr*pi;
nts=round(ts/dt);
dx=1;
parfor i=1:ntr
    % receiver location
    rr = R;
    thetar = 0.5*pi;
    phir =phiar(i);
    xr=rr*sin(thetar)*cos(phir);
    yr=rr*sin(thetar)*sin(phir);
    zr=rr*cos(thetar);
    
    [urdum,usynt,usynp,dtm,ntm] =syn_modes_func_cart(...
        modesS,modesT,xs,ys,zs,rr,thetar,phir,...
        force_convert(fsrc,thetas,phis),vp,vs,Q,dt,nt);
    urdumsum=urdum;
    uzmode_g(:,i)=urdumsum.*exp(-(0:nt-1)*dt*4/T);
    umoderw=ifft(urdumsum.*exp(-(0:nt-1)*dt*4/T))*dtm*ntm;
    umodet=real(fft(umoderw.*wf))*dw/(2*pi);
    uzmode(:,i)=umodet;
end
%%
nfold=20;
figure
wiggle(time-ts,(1:36)*5,uzmode.',1,'p',nfold)
xlim([0 40])
%%
save('modes_syn_data20km_sgf_dep10kmQ50.mat','time','uzmode','uzmode_g');
%save('modes_syn_data20km_sgf3.mat','time','uzmode');
%save('modes_syn_data20km_sgfz_0p3pi.mat','time','uzmode');
%save('modes_syn_data20km_sgfz_dep4km.mat','time','uzmode','uzmode_g')
%%
f0=0.5;
ts =1.5/f0; 
 [w,t] = ricker(f0,ts,dt,nt); 
 wf = ifft(w)*dt*nt;
for i=1:ntr
    umoderw=ifft(uzmode_g(:,i).'/dx)*dt*nt;
    umodet=real(fft(umoderw.*wf))*dw/(2*pi);
    uzmode(:,i)=umodet;
end