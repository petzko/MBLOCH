function erg=gain(path,ch,d,fmin,fmax,sw,ni1,nf1)
%sw==1 plot oscillator strength*relative inversion, sw==2 plot oscillator strength*upper level occupation, sw==0 plot gain curve
%d is length of one period in Angstrom
%fmin and fmax are the minimum and maximum considered frequency in Hz
n0=3.8; meff=0.067; %meff only needed if sw==1, n0 only needed for the other cases
eps0=8.8541878176e-12; e0=1.60217646e-19; hbar=1.05457168e-34; c0=3e8; me=9.10938188e-31;
global p; global n;
[dummy,ve]=sort(ch); m=length(ve); %index matrix ind
if(~exist('sw')) sw=0; end; %if(sw~=1) sw=0; end;
if(~exist('ni1')) ni1=0; else ni1h=2*m+mod(ni1-1,m)+1; nf1=ni1h+nf1-ni1; ni1=ni1h; end;
fid=fopen([path,'/gain.dat']);
if(fid~=-1)
  fclose(fid);
  A=load([path,'/gain.dat']);
  E=A(:,1); h=(A(:,(2*m+2):(4*m+1))+A(:,(4*m+2):(6*m+1)))/2; n=h(:,1:2:(2*m)); p=h(:,2:2:(2*m)); %3d density of particles in a state n in energy cell around E(i) is given by A(i,n+1)*|psi_nsu|^2
    
  psief=load(strcat(path,'\psi.dat'));
  eef=load(strcat(path,'\e_bound'));
  nz=length(psief(:,2))/(4*m);
  z=psief(1:nz,1);
  ez=zeros(4*m,4*m); wif=zeros(4*m,4*m);
  for ni=(2*m+1):(3*m)
      for nf=1:(4*m)
          wif(ni,nf)=max(0,(eef(ni,2)-eef(nf,2))*e0/hbar);
          if(wif(ni,nf)>0)
              psieu=psief(((ni-1)*nz+1):(ni*nz),2); psieu=psieu/sqrt(trapz(z,psieu.^2));
              psiel=psief(((nf-1)*nz+1):(nf*nz),2); psiel=psiel/sqrt(trapz(z,psiel.^2));
              
              %calculate the dipole element! 
              zs=trapz(z,z.*abs(psiel.*psieu)); 
              if(zs ~=0) zs = zs/trapz(z,abs(psiel.*psieu)); end;
              ez(ni,nf)=e0*trapz(z-zs,psiel.*(z-zs).*psieu)*1e-10;
              
              
              %ez=trapz(interp(z,4),interp(psiel,4).*interp(z,4).*interp(psieu,4));
          end;
      end;
  end;
  %wmin=min(min(wif.*(1+0./wif))); wmax=max(max(wif));
  k=0;
  for w=2*pi*(fmin+(fmax-fmin)*(1:100)/100)
      k=k+1; f(k)=w/2/pi; g(k)=0;
      for ni=(2*m+1):(3*m)
          for nf=1:(4*m)
              if((ni1==0)||((ni==ni1)&&(nf==nf1)))
              dab=ez(ni,nf); wab=wif(ni,nf);
              if(wab>0)
                  for nE=(1:length(E))
                      nih=mod(ni-1,m)+1; nfh=mod(nf-1,m)+1; Delta=(n(nE,nih)-n(nE,nfh))/(d*1e-10);
                      if(n(nE,nih)*p(nE,nih)==0) if(nE>1)gi=gi_old; else gi=0; end; else gi=(p(nE,nih)/n(nE,nih)); end;
                      if(n(nE,nfh)*p(nE,nfh)==0) if(nE>1)gf=gf_old; else gf=0; end; else gf=(p(nE,nfh)/n(nE,nfh)); end;
                      gi_old=gi; gf_old=gf;
                      gamma=(gi+gf)/2; %assumption that in stationary case, in- and out-scattering rates are equal (if we only consider intersubband scattering, this does not hold for each individual transition! However, the numerical example GaAs3.4_15test yielded practically identical results for the spectral gain obtained based on in- and outscattering)
                      if(gamma~=0) g(k)=g(k)+w/c0/n0*dab^2/gamma/eps0/hbar*Delta/(1+(w-wab)^2/gamma^2); end; %power gain
                  end;
              end;
              end;
          end;
      end;
  end;
  %if(~sw) plot(f/1e12,g/100); erg=[f;g]; [maxg,imax]=max(g); disp(['f=',num2str(f(imax)/1e12),' THz, FWHM=',num2str(fw(g)*(f(2)-f(1))/1e12),' THz, max=',num2str(maxg/100),' 1/cm']); end;
  if(~sw) 
      plot(f/1e12,g/100); 
      erg=[f;g]; 
      [maxg,imax]=max(g); 
      disp(['f=',num2str(f(imax)/1e12),' THz, max=',num2str(maxg/100),' 1/cm']); 
  end;

  k=0; occ=sum(n); occ=occ/sum(occ); 
  for ni=(2*m+1):(3*m)
      for nf=1:(4*m)
          if((ni1==0)||((ni==ni1)&&(nf==nf1)))
          wab=wif(ni,nf);
          nih=mod(ni-1,m)+1; nfh=mod(nf-1,m)+1;
          if(wab>0)
              k=k+1; fab(k)=wab/2/pi;
              fosc(k)=2*(ez(ni,nf)/e0)^2*me*meff*wab/hbar;
              if(sw==1)fd(k)=fosc(k)*(occ(nih)-occ(nfh)); end;
              if(sw==2)fd(k)=fosc(k)*occ(nih); end;
              stat(k,1)=nih; stat(k,2)=nfh;
          end;
          end;
      end;
  end;
  if(sw==1) plot(fab,fd,'rx','markersize',12); sp=max(fd)-min(fd); if(sp>0) axis([fmin fmax min(fd)-0.1*sp max(fd)+0.1*sp]); end; xlabel('Frequency [Hz]'); ylabel('\Delta_p*f_{osc}'); erg=[fab;fd];
      [dummy,ind]=sort(fd); ind=fliplr(ind); disp('f_osc*Delta   f [THz]     n_i       n_f'); mi=min(length(ind),4); disp([fd(ind(1:mi));fab(ind(1:mi))/1e12;ch(stat(ind(1:mi),:))']');
  end;
  if(sw==2) plot(fab,fd,'rx','markersize',12); sp=max(fd)-min(fd); if(sp>0) axis([fmin fmax min(fd)-0.1*sp max(fd)+0.1*sp]); end; xlabel('Frequency [Hz]'); ylabel('p_i*f_{osc}'); erg=[fab;fd];
      [dummy,ind]=sort(fd); ind=fliplr(ind); disp('f_osc*p_i   f_osc     f [THz]     p_i       p_f'); mi=min(length(ind),4); disp([fd(ind(1:mi));fosc(ind(1:mi));fab(ind(1:mi))/1e12;ch(stat(ind(1:mi),:))']');
  end;

else erg=ch*NaN; end