function [erg,t,I]=currPETZ(path,ave,tg,plotonoff)
curr=load([path,'\dens_corrente.dat']);
n=length(curr(:,1));
t=curr(:,1); I=smooth(curr(:,2),ave); p=(1:(n-ave))+ave/2;
if(plotonoff)
plot(t(p)*1e12,I(p));
end
xlabel('t/ps'); ylabel('J/(A/m^2)');
tg=tg/1e15; erg=sum((curr(:,1)>tg).*curr(:,2))/sum((curr(:,1)>tg))/1e7