function [E,dis]=plotdistr(ave,ch)
%e.g., ch=[3 4 1 5 2] means that the level with the second-lowest energy
%should be labeled '4' (and obtain the respective color), which might indicate the lower laser
%level, and so on
%path='C:\jirausch\tum\sige\ge\twobandge2_110';
%erg=fde(tg_in_fs,path);
%plot(ave(:,1),ave(:,20:23)); legend('1','2','3','4','5','6','7','8'); xlabel('E/eV');
%hold off;
global T; global T2; global occ; global dis0; T=[]; T2=[]; occ=[]; dis0=[];
if ischar(ave) path1=ave; ave=load([path1,'\fe.dat']); end;
e0=1.60217646e-19; kB=1.3806503e-23;
plot(0,0); hold on;
B =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0;0.5,0.5,0.5;0,0,0;0.1,0.1,0.1;0.9,0.9,0.9];
m=length(ch);
A=B(1:m,:);
set(gca,'ColorOrder',A);
[dummy,ind] = sort(ch); %index matrix ind
E=ave(1:(length(ave(:,1))),1); dis=ave(1:(length(ave(:,1))),1+m+ind); %sort curves
%E=E(1:(end-1)); dis=dis(1:(end-1),:); %sort curves
plot(E,dis(:,:),'linewidth',2);
legendinfo = {'dummy'}; 
for(i=1:length(ch))
    legendinfo{i+1} = num2str(ch(i));
end
legend(legendinfo);

%Epop=0.036;
%plot([E-0.0036 E+0.0027 E-0.0174 E-0.0131 E],dis,'linewidth',2);
%plot([E-0.0036 E+0.0027 E-0.0174 E-0.0131 E]+0.0571,dis,'linewidth',2);
xlabel('E/eV'); ylabel('f(E)');
T=sum((E*ones(1,m)).*dis)*e0./sum(dis)/kB
%for n=1:m [a(n),resnorm(n)]=lsqcurvefit(@fE,100,E-E(1),dis(:,n)/dis(1,n)); end;
%for n=1:m [x,resnorm(n)]=lsqcurvefit(@fE1,[dis(1,n) e0/kB/T(n)],E,dis(:,n)); A0(n)=x(1); T2(n)=e0/kB/x(2); end;
%for n=1:m plot(E,dis(1,n)*fE(a(n),E-E(1)),':'); end;
%for n=1:m plot(E,fE1([A0(n) e0/kB/T2(n)],E),':'); end;
hold off;
T
T2
%erg=300./(0.026*a)
% resnorm
%sum(dis(:,1:m))/sum(sum(dis(:,1:m)))
dis1=ave(1:(end-0),1+m+ind); occ1=sum(dis1(:,1:m))/sum(sum(dis1(:,1:m)))
dis2=ave(1:(end-0),1+2*m+ind); occ2=sum(dis2(:,1:m))/sum(sum(dis2(:,1:m)))
dis0=(dis1+dis2)/2; occ=sum(dis0(:,1:m))/sum(sum(dis0(:,1:m)))
share=sum(dis(end,:))/sum(sum(dis(:,:)))
if exist(path1) save([path1,'\occ'],'occ','-ASCII'); end; %occ nochmal checken; richtige Zuordnung in ps Programm (kein ch); Datei schreiben nur bei Bedarf
%erg=sum(dis(:,1:m)+dis2(:,1:m))/sum(sum(dis(:,1:m)+dis2(:,1:m)));