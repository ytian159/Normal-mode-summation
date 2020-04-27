function [usynr,usynt,usynp,dt,npts] =syn_modes_func3(modesS,modesT,rs,thetas,phis,rr,thetar,phir,f0,vp,vs,Q,dt,npts)
nl=length(modesS);
% dt=40;
% npts=1000;
t=(0:npts-1)*dt;
usynr=0;
usynt=0;
usynp=0;
for k=1:nl
    ao = modesS(k).L;
    l=ao;
%     if l==1
%         continue;
%     end
    mfreq = modesS(k).freq;
    nmls=length(mfreq);
    mfreqt=modesT(k).freq;
    nmlt=length(modesT(k).freq);
    rat1S=modesS(k).rat1;
    rat2S=modesS(k).rat2;
    nmS=modesS(k).fac;
    nmT=modesT(k).fac;
    for j=1: nmls
        f = mfreq(j);
        w = 2*pi*f;
        ka = w/vp;
        kb = w/vs;
        
        [Us,Vs]=RadSphMod2(l,ka,kb,rs,rat1S(j),rat2S(j));
        [Ur,Vr]=RadSphMod2(l,ka,kb,rr,rat1S(j),rat2S(j));
        Usn=Us*nmS(j);
        Vsn=Vs*nmS(j);
        Urn=Ur*nmS(j);
        Vrn=Vr*nmS(j);
        %tfac=(1-exp(-w*t/(2*Q)).*cos(w*t))/w^2;
        tfac=sin(w*t).*exp(-w*t/(2*Q))/w;
        subur=0;
        subut=0;
        subup=0;
        
        for m=-ao:ao
            [Rrr,Rrt,Rrp]=Rlm_func(l,m,thetar,phir);
            [Srr,Srt,Srp]=Slm_func(l,m,thetar,phir);
            %[Trr,Trt,Trp]=Tlm_func(l,m,thetar,phir);
            
            [Rsr,Rst,Rsp]=Rlm_func(l,m,thetas,phis);
            [Ssr,Sst,Ssp]=Slm_func(l,m,thetas,phis);
            %[Tsr,Tst,Tsp]=Tlm_func(l,m,thetas,phis);
            
            usr=Usn*Rsr+Vsn*Ssr;
            ust=Usn*Rst+Vsn*Sst;
            usp=Usn*Rsp+Vsn*Ssp;
            
            urr=Urn*Rrr+Vrn*Srr;
            urt=Urn*Rrt+Vrn*Srt;
            urp=Urn*Rrp+Vrn*Srp;
            
            
            Amp=conj(usr*f0(1)+ust*f0(2)+usp*f0(3));
            
            subur=subur+Amp*urr;
            subut=subut+Amp*urt;
            subup=subup+Amp*urp;
            
            
        end
        
        usynr=usynr+subur*tfac;
        usynt=usynt+subut*tfac;
        usynp=usynp+subup*tfac;
        
    end
    if l==0
        continue;
    end
    for j=1: nmlt
        f = mfreqt(j);
        
        w = 2*pi*f;
        %ka = w/vp;
        kb = w/vs;
        Ws=RadTorMod(l,kb,rs);
        Wr=RadTorMod(l,kb,rr);
        Wsn=Ws*nmT(j);
        Wrn=Wr*nmT(j);
        %tfac=(1-exp(-w*t/(2*Q)).*cos(w*t))/w^2;
        tfac=sin(w*t).*exp(-w*t/(2*Q))/w;
        subur=0;
        subut=0;
        subup=0;
        for m=-ao:ao
            [Trr,Trt,Trp]=Tlm_func(l,m,thetar,phir);
            [Tsr,Tst,Tsp]=Tlm_func(l,m,thetas,phis);
            
            usr=Wsn*Tsr;
            ust=Wsn*Tst;
            usp=Wsn*Tsp;
            
            urr=Wrn*Trr;
            urt=Wrn*Trt;
            urp=Wrn*Trp;
            
            Amp=conj(usr*f0(1)+ust*f0(2)+usp*f0(3));
            
            subur=subur+Amp*urr;
            subut=subut+Amp*urt;
            subup=subup+Amp*urp;
        end
        
        usynr=usynr+subur*tfac;
        usynt=usynt+subut*tfac;
        usynp=usynp+subup*tfac;
        
    end
    
    
end

end