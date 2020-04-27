%%
clearvars;
set(0,'defaultaxesfontsize',16); % control the default fontsize of the plots
set(0,'defaulttextfontsize',16);

vp = 6e3;
vs = 3e3;

rho = 3000;
mu = rho* vs^2;
lamda = rho*vp^2 - 2*mu;
R= 20e3;
nr = 2e3;
ri = linspace(1, R, nr);
%nr=length(ri);
dr=ri(2)-ri(1);
angord =[0 1:60];
freq = [0.01:0.01/10:1.5];
nf=length(freq);
nl=numel(angord);
equ_num(1:nl)=0;
equT(1:nl)=0;
% x1=zeros(nl,1);
% x2=zeros(nl,1);

%%
tic
for imod = 1: numel(angord)
    l = angord(imod)
    equ_num=zeros(nf,1);
    
    for j=1: nf
        f = freq(j);
        w = 2*pi*f;
        ka = w/vp;
        kb = w/vs;
        CRcc=ka.^(-1).*kb.^(-1).*R.^(-2).*((-2).*kb.*l.*mu.*SphericalBesselJ(l, ...
            ka.*R)+2.*kb.*l.^2.*mu.*SphericalBesselJ(l,ka.*R)+(-1).*ka.^2.* ...
            kb.*lamda.*R.^2.*SphericalBesselJ(l,ka.*R)+(-2).*ka.^2.*kb.*mu.* ...
            R.^2.*SphericalBesselJ(l,ka.*R)+4.*ka.*kb.*mu.*R.* ...
            SphericalBesselJ(1+l,ka.*R));
        CRbb=ka.^(-1).*kb.^(-1).*R.^(-2).*((-2).*ka.*l.*mu.*SphericalBesselJ(l, ...
            kb.*R)+2.*ka.*l.^3.*mu.*SphericalBesselJ(l,kb.*R)+(-2).*ka.*kb.* ...
            l.*mu.*R.*SphericalBesselJ(1+l,kb.*R)+(-2).*ka.*kb.*l.^2.*mu.*R.* ...
            SphericalBesselJ(1+l,kb.*R));
        CScc=ka.^(-1).*kb.^(-1).*mu.*R.^(-2).*((-2).*kb.*SphericalBesselJ(l, ...
            ka.*R)+2.*kb.*l.*SphericalBesselJ(l,ka.*R)+(-2).*ka.*kb.*R.* ...
            SphericalBesselJ(1+l,ka.*R));
        CSbb=ka.^(-1).*kb.^(-1).*mu.*R.^(-2).*((-2).*ka.*SphericalBesselJ(l, ...
            kb.*R)+2.*ka.*l.^2.*SphericalBesselJ(l,kb.*R)+(-1).*ka.*kb.^2.* ...
            R.^2.*SphericalBesselJ(l,kb.*R)+2.*ka.*kb.*R.*SphericalBesselJ(1+ ...
            l,kb.*R));
        % x1(j)=-CScc/CSbb;
        %         x2(j)=-CRcc/CRbb;
        equ_num(j)=CRcc*CSbb-CRbb*CScc;
        %equ_num(j)=-CScc/CSbb+CRcc/CRbb;
        equT(j)=mu.*R.^(-1).*(((-1)+l).*SphericalBesselJ(l,kb.*R)+(-1).*kb.*R.* ...
            SphericalBesselJ(1+l,kb.*R));
        if l==0
            %             tmp1 = (ka*R*(lamda+2*mu)*SphericalBesselJ(0,ka*R)+(3*lamda-2*mu)*SphericalBesselJ...
            %                 (1,ka*R)-ka*R*(lamda+2*mu)*SphericalBesselJ(2,ka*R));
            
            %             tmp1=(-1).*ka.^(-1).*kb.^(-1).*mu.*R.^(-2).*((-2).*kb.* ...
            %                 SphericalBesselJ(0,ka.*R)+(-2).*ka.*kb.*R.*SphericalBesselJ(1,ka.* ...
            %                 R));
            tmp1=ka.^(-1).*kb.^(-1).*R.^(-2).*((-1).*ka.^2.*kb.*(lamda+2.*mu).* ...
                R.^2.*SphericalBesselJ(0,ka.*R)+4.*ka.*kb.*mu.*R.* ...
                SphericalBesselJ(1,ka.*R));
            equ_num(j) = tmp1;
        end
    end
    % find zeros (mode frequecies)
    k = 0;
    xzero=[];
    if nf < 2
        continue;
    end
    
    for j =1: nf-1
        if equ_num(j)*equ_num(j+1) <= 0
            if equ_num(j)==equ_num(j+1)
                continue;
            end
            xtmp=interp1([equ_num(j) equ_num(j+1)],[freq(j) freq(j+1)],0);
            if isnan(xtmp)
                continue;
            end
            if xtmp<l*3e-3
                continue;
            end
            if ~isempty (xtmp)
                k = k + 1;
                
                %                 rat1tmp=interp1([freq(j) freq(j+1)],[x1(j) x1(j+1)],xtmp)
                %                 rat2tmp=interp1([freq(j) freq(j+1)],[x2(j) x2(j+1)],xtmp)
                xzero(k)=xtmp;
            end
            
        end
    end
    if l==1
        xzero(1)=[];
    end
    %%
    xtzero=[];
    k2=0;
    
    for j =1: nf-1
        if equT(j)*equT(j+1) <= 0
            if equT(j)==equT(j+1)
                continue;
            end
            xtmpT=interp1([equT(j) equT(j+1)],[freq(j) freq(j+1)],0);
            if isnan(xtmpT)
                continue;
            end
            if xtmpT<l*3e-3
                continue;
            end
            if ~isempty (xtmpT)
                k2 = k2 + 1;
                if l==1 && k2==1
                    xtzero(k2)=0;
                    continue
                end
                xtzero(k2)=xtmpT;
            end
            
        end
    end
    if l==0
        xtzero=[];
    end
    if l==1
        xtzero(1)=[];
    end
    modesS(imod).L = l;
    modesS(imod).freq = xzero;
    
    modesT(imod).L = l;
    modesT(imod).freq = xtzero;
end
toc
%%
figure;
hold on;
for k=1: numel(angord)
    if isempty(modesS(k).freq)
        continue;
    end
    plot( modesS(k).L, modesS(k).freq*1e3,'ko','markerfacecolor','r');
end
figure;
hold on;
for k=1: numel(angord)
    if isempty(modesT(k).freq)
        continue;
    end
    plot( modesT(k).L, modesT(k).freq*1e3,'ko','markerfacecolor','r');
end
pause(0.1)
%%

for k=1:nl
    ao = modesS(k).L;
    l=ao;
    mfreq = modesS(k).freq;
    nml=length(mfreq);
    nmSfac=zeros(nml,1);
    nmTfac=zeros(nml,1);
    x1=zeros(nml,1);
    x2=zeros(nml,1);
    mfreqt=modesT(k).freq;
    nmlt=length(modesT(k).freq);
    for j=1: nml
        f = mfreq(j);
        w = 2*pi*f;
        ka = w/vp;
        kb = w/vs;

        Sa=SphericalBesselJ(l+1,ka*R)/SphericalBesselJ(l,ka*R);
        Sb=SphericalBesselJ(l+1,kb*R)/SphericalBesselJ(l,kb*R);
        x1(j)=-kb*(2*l*(l-1)-(ka*R*kb/ka)^2+4*ka*R*Sa)*SphericalBesselJ(l,ka*R)/(2*ka*l*(l+1)*(l-1-kb*R*Sb)*SphericalBesselJ(l,kb*R));
        x2(j)=2*kb*(1-l+ka*R*Sa)*SphericalBesselJ(l,ka*R)/(ka*SphericalBesselJ(l,kb*R)*(-2+2*l^2-(kb*R)^2+2*kb*R*Sb));
        if l==0
            x1(j)=x2(j);
        end
        [U,V]=RadSphMod2(l,ka,kb,ri,x1(j),x2(j));
        
        %[U,V]=radialmodes ( l,vp, vs,w,x1(j), ri);
        %         figure
        %         plot(Ur);
        
        %         U(:,j)=Ur;
        %         V(:,j)=Vr;
        %         W(:,j)=Wr;
        nmSmod=(U.^2+V.^2).*ri.^2*rho*dr;
        nmSfac(j)=sqrt(1./sum(nmSmod));
        
    end
    
    for j=1: nmlt
        f = mfreqt(j);
        
        w = 2*pi*f;
        %ka = w/vp;
        kb = w/vs;
        [Wr]=RadTorMod(l,kb,ri);
        nmTfac(j)=1/sqrt(sum(Wr.^2.*ri.^2*dr*rho));
    end
    %x3=(x1-x2)./x1
    %x1
    
    modesS(k).rat1=x1;
    modesS(k).rat2=x2;
    modesS(k).fac=nmSfac;
    if ao==0
        nmTfac=0;
    end
    modesT(k).fac =nmTfac;
    %dum=max((x1-x2)/x1)
end
%%
%save('modes_data_20km_1p5hz.mat','modesS','modesT');