function erg=scan(path_,ni,nf,C,E,sw)
%for investigation of a single structure set sw=3 and enter complete path of this structure
%T38=scan('f:/jirauschek/gaas3.8_',5,3,C38,E38,0);

%new MC (with impurity scattering)
%1: e-e scattering: 2: optical absorption; 3: optical emission; 4+5 reserved for other absorption/emission processes between same type of valley;
%6: acoustic phonons; 7: impurity; 8: interface roughness; 9-11: scattering between different types of valleys

%old MC files (without impurity scattering)
%1: polar abs; 2: polar em.; 3: acoustic; 4: non-polar intervalley absorpt. (diff. type of valleys, e.g. Gamma-L);
%5: non-polar intervalley em. (diff. type of valleys, e.g. Gamma-L); 6: carrier-carrier;
%7: Interband intrasubband (i.e., both carriers remaining in their initial subbands) carrier-carrier scattering
%8: non-polar intervalley absorpt. (same type of valleys in cond. band, e.g. L-L)
%9: non-polar intervalley em. (same type of valleys in cond. band, e.g. L-L)
%10: non-polar intravalley absorption (for holes)
%11: non-polar intravalley emission (for holes)
m=length(C(1,:));
global dE; global rf; global rb; global dEh;
global sr_op; global sr_ac; global sr_ee; global sr_opem; global sr_opabs; global sr_imp; global sr_if; global sr_las; global sr_sp; global sr_al;
global sf_op; global sf_ac; global sf_ee; global sf_opem; global sf_opabs; global sf_imp; global sf_if; global sf_las; global sf_sp; global sf_al; global sf;
global sb_op; global sb_ac; global sb_ee; global sb_opem; global sb_opabs; global sb_imp; global sb_if; global sb_las; global sb_sp; global sb_al; global sb;

dE=0; rf=0; rb=0; dEh=0; sr_op=0; sr_ac=0; sr_ee=0; sr_imp=0; sr_if=0; sr_las=0; sr_sp=0; sr_al=0; sr_opem=0; sr_opabs=0; sf_op=0; sf_ac=0;
sf_ee=0; sf_imp=0; sf_if=0; sf_las=0; sf_sp=0; sf_al=0; sf_opem=0; sf_opabs=0; sf=0; sb_op=0; sb_ac=0; sb_ee=0; sb_imp=0; sb_if=0; sb_las=0; sb_sp=0; sb_al=0; sb_opem=0; sb_opabs=0; sb=0;

dEh=0;
if(sw==3) nv=1; else nv=1:length(E); end;
for n=nv
   if(sw==3) path=path_; else path=[path_,num2str(n)]; end;
   A_=load([path,'/srates.dat']); mat_=load([path,'/srateprep.dat']); E_=load([path,'/e_bound']);
   sr(n)=0; sr_op(n)=0; sr_ac(n)=0; sr_ee(n)=0; sr_imp(n)=0; sr_if(n)=0; sr_opem(n)=0; sr_opabs(n)=0; sr_las(n)=0; sr_sp(n)=0; sr_al(n)=0;
   sf(n)=0; sf_op(n)=0; sf_ac(n)=0; sf_ee(n)=0; sf_imp(n)=0; sf_if(n)=0; sf_opem(n)=0; sf_opabs(n)=0; sf_las(n)=0; sf_sp(n)=0; sf_al(n)=0;
   sb(n)=0; sb_op(n)=0; sb_ac(n)=0; sb_ee(n)=0; sb_imp(n)=0; sb_if(n)=0; sb_opem(n)=0; sb_opabs(n)=0; sb_las(n)=0; sb_sp(n)=0; sb_al(n)=0;
   for k=1:length(nf)
       if(n>nv(1))
        srates2(A_,mat_,E_,ni,nf(k),C(n,:),-1); dE1(1)=dE; srates2(A_,mat_,E_,ni,nf(k)-m,C(n,:),-1); dE1(2)=dE; srates2(A_,mat_,E_,ni,nf(k)+m,C(n,:),-1); dE1(3)=dE;
        [dummy,ind] = min(abs(dE1-dEh(n-1,k)));
       else ind=1; end;
       ni1=ni; nf1=nf(k); %if(ind==2) nf1=nf(k)-m; elseif(ind==3) nf1=nf(k)+m; end;
       dn(k)=ni1-nf1;
       if((k==1)|(min(abs(dn(1:(k-1))-dn(k)))>0)) %avoid that a level is counted more than once
        sr(n)=sr(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),0); dEh(n,k)=dE; sf(n)=sf(n)+rf; sb(n)=sb(n)+rb;
        if(sw~=2)
            sr_opem(n)=sr_opem(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),3); sf_opem(n)=sf_opem(n)+rf; sb_opem(n)=sb_opem(n)+rb;
            %sr_opem(n)=sr_opem(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),2); sf_opem(n)=sf_opem(n)+rf; sb_opem(n)=sb_opem(n)+rb;
            %sr_opem(n)=sr_opem(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),5); sf_opem(n)=sf_opem(n)+rf; sb_opem(n)=sb_opem(n)+rb;
            %sr_opem(n)=sr_opem(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),9); sf_opem(n)=sf_opem(n)+rf; sb_opem(n)=sb_opem(n)+rb;
            %sr_opem(n)=sr_opem(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),11); sf_opem(n)=sf_opem(n)+rf; sb_opem(n)=sb_opem(n)+rb;
            
            sr_opabs(n)=sr_opabs(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),2); sf_opabs(n)=sf_opabs(n)+rf; sb_opabs(n)=sb_opabs(n)+rb;
            %sr_opabs(n)=sr_opabs(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),1); sf_opabs(n)=sf_opabs(n)+rf; sb_opabs(n)=sb_opabs(n)+rb;
            %sr_opabs(n)=sr_opabs(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),4); sf_opabs(n)=sf_opabs(n)+rf; sb_opabs(n)=sb_opabs(n)+rb;
            %sr_opabs(n)=sr_opabs(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),8); sf_opabs(n)=sf_opabs(n)+rf; sb_opabs(n)=sb_opabs(n)+rb;
            %sr_opabs(n)=sr_opabs(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),10); sf_opabs(n)=sf_opabs(n)+rf; sb_opabs(n)=sb_opabs(n)+rb;
            
            sr_op(n)=sr_opem(n)+sr_opabs(n); sf_op(n)=sf_opem(n)+sf_opabs(n); sb_op(n)=sb_opem(n)+sb_opabs(n);
           
            sr_ac(n)=sr_ac(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),6); sf_ac(n)=sf_ac(n)+rf; sb_ac(n)=sb_ac(n)+rb;
            %sr_ac(n)=sr_ac(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),3); sf_ac(n)=sf_ac(n)+rf; sb_ac(n)=sb_ac(n)+rb;
            sr_ee(n)=sr_ee(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),1); sf_ee(n)=sf_ee(n)+rf; sb_ee(n)=sb_ee(n)+rb;
            %sr_ee(n)=sr_ee(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),6); sf_ee(n)=sf_ee(n)+rf; sb_ee(n)=sb_ee(n)+rb;
            %sr_ee(n)=sr_ee(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),7); sf_ee(n)=sf_ee(n)+rf; sb_ee(n)=sb_ee(n)+rb;
            sr_imp(n)=sr_imp(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),7); sf_imp(n)=sf_imp(n)+rf; sb_imp(n)=sb_imp(n)+rb;
            sr_if(n)=sr_if(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),8); sf_if(n)=sf_if(n)+rf; sb_if(n)=sb_if(n)+rb;
            sr_las(n)=sr_las(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),9); sf_las(n)=sf_las(n)+rf; sb_las(n)=sb_las(n)+rb;
            sr_sp(n)=sr_sp(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),4); sf_sp(n)=sf_sp(n)+rf; sb_sp(n)=sb_sp(n)+rb;
            sr_al(n)=sr_al(n)+srates2(A_,mat_,E_,ni1,nf1,C(n,:),5); sf_al(n)=sf_al(n)+rf; sb_al(n)=sb_al(n)+rb;
       end;
       end;
   end;
end;
if(sw==1) subplot(2,1,1); set(gca,'NextPlot','replacechildren'); A =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0]; [dummy,nfind] = sort(nf);
    B=A(mod(nf-1,m)+1,:); set(gca,'ColorOrder',B); plot(E,dEh); subplot(2,1,2); end;
if((sw==1)|(sw==0))
    plot(E,sr,'b',E,sr_op,'b:',E,sr_ac,'b-.',E,sr_ee,'b--',sr_imp,'b-x',sr_if,'b-+'); xlabel('E/(kV/cm)'); ylabel('rate/s^-^1'); end;
if(sw==3) disp([sr_op sr_ac sr_ee sr_imp sr_if sr_las sr_sp]); end;
%sr_las
%sf_sp
%legend('total','optical phonons','acoustic phonons','e-e scattering','impurity');
erg=sr;
%erg=round([sr_op sr_ac sr_ee sr_imp sr_if]/1e9);