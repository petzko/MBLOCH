function res=srates2(A,mat,E,ni,nf,ch,itype)
%new MC (with impurity scattering)
%1: e-e scattering: 2: optical absorption; 3: optical emission; 4+5 reserved for other absorption/emission processes between same type of valley;
%6: acoustic phonons; 7: impurity; 8: reserved for interface roughness or similar; 9-11: scattering between different types of valleys
%here: 0: total scattering; -1: only dE

%old MC files (without impurity scattering)
%1: polar abs; 2: polar em.; 3: acoustic; 4: intervalley absorpt. (diff. type of valleys, e.g. Gamma-L);
%5: intervalley em. (diff. type of valleys, e.g. Gamma-L);
%6: carrier-carrier; 7: intervalley absorpt. (same type of valleys, e.g. L-L) or for holes transverse optical absorpt;
%8: intervalley em. (same type of valleys, e.g. L-L) or for holes transverse optical emiss)

nsub=sqrt(length(A)/11); nperiod=nsub/4; Bf=reshape(A,nsub,nsub,11);
C=mat((length(mat)-length(E)+1):length(mat));
global rf; global rb;

%to obtain the effective scattering rate
for ni1=(nperiod+1):(3*nperiod) for nf1=1:nsub
   nf1eq=nf1; if(nf1eq>3*nperiod) nf1eq=nf1eq-2*nperiod; elseif(nf1eq<=nperiod) nf1eq=nf1eq+2*nperiod; end; ni2=nf1eq; %equivalent levels
   nf2=ni2+ni1-nf1; if(nf2>nsub) ni2=ni2-nperiod; nf2=nf2-nperiod; elseif(nf2<=0) ni2=ni2+nperiod; nf2=nf2+nperiod; end;
   Bb(nf1,ni1,:)=C(nf1eq)/C(ni1)*Bf(nf2,ni2,:);
   Beff(nf1,ni1,:)=Bf(nf1,ni1,:)-Bb(nf1,ni1,:);
end; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dummy,ind] = sort(ch); %index matrix ind
m=length(ch);
%ind(1:2)=ind(1:2)-m*(ind(1:2)>ind(3)); ind(4:5)=ind(4:5)+m*(ind(4:5)<ind(3)); if(min(ind)<1) ind=ind+m; end;  %assumption that within one period, 1,2 are underneath 3 and 4,5 are above 3
if((ni>-100)&(nf>-100))
    nf=m+ind(mod(nf-1,m)+1)+m*floor((nf-1)/m)-m*floor((ni-1)/m); ni=m+ind(mod(ni-1,m)+1);
    if(ni>nf) dn=m*ceil(ni/m-3); ni=ni-dn; nf=nf-dn; end; %shift ni to 3rd period
    if(ni<=nf) dn=m*ceil(ni/m-2); ni=ni-dn; nf=nf-dn; end; %shift ni to 2nd period
    ni_=ni; nf_=nf; %average
    if((ni<=2*nperiod)&(nf<=3*nperiod)) ni_=ni+nperiod; nf_=nf+nperiod;
    elseif((ni>2*nperiod)&(nf>nperiod)) ni_=ni-nperiod; nf_=nf-nperiod;
    end;
    global dE; global fif;
    nie=mod(ni-1,nperiod)+1; nfe=mod(nf-1,nperiod)+1; dE=(E(nfe,2)-E(nie,2))/(E(m+1,2)-E(1,2))+floor((nf-1)/nperiod)-floor((ni-1)/nperiod); fif=(E(nie,2)-E(nfe,2))*2.4180e14;
    if((nf<1)|(nf>(4*nperiod))) res=0; rf=0; rb=0; else
        if(itype==0) res=0.5*sum(Beff(nf,ni,1:11)+Beff(nf_,ni_,1:11)); rf=0.5*sum(Bf(nf,ni,1:11)+Bf(nf_,ni_,1:11)); rb=0.5*sum(Bb(nf,ni,1:11)+Bb(nf_,ni_,1:11));
        elseif(itype<0);
        else res=0.5*(Beff(nf,ni,itype)+Beff(nf_,ni_,itype)); rf=0.5*(Bf(nf,ni,itype)+Bf(nf_,ni_,itype)); rb=0.5*(Bb(nf,ni,itype)+Bb(nf_,ni_,itype)); end;
    end;
elseif(ni>-100)
    ni=m+ind(mod(ni-1,m)+1)+m*floor((ni-1)/m); nf=[(1:(ni-1)),((ni+1):nsub)];
    if(itype==0) res=sum(sum(Beff(nf,ni,1:11))); rf=sum(sum(Bf(nf,ni,1:11))); rb=sum(sum(Bb(nf,ni,1:11)));
    else res=sum(Beff(nf,ni,itype)); rf=sum(Bf(nf,ni,itype)); rb=sum(Bb(nf,ni,itype)); end;
end;
if(itype>=0) if((~(res>0))&(~(res<=0))) res=0; rf=0; rb=0; end; end; %replace NaN by 0