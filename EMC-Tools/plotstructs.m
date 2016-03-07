%**************************************************************************
%
% This programm plots the structure of the conduction band edge and the
% squared absolute values of the wavefunctions for the calculated bound
% state energy levels
%
% FUNCTION []=
%         PLOTSTRUCT(DIRECTORY)
%
%DIRECTORY....................the directory in the home directory in which
%.............................the structure files are stored 
%.............................(e.g. 'Project\QCL')
%
%**************************************************************************

function [] =...
	plotstructs(k,dir,ch)
%e.g., ch=[3 4 1 5 2] means that the level with the second-lowest energy
%should be colored with color '4', which might indicate the lower laser
%level, and so on
global mpt;
fa=3.0;
ppot=load (strcat(dir,'\ppot',num2str(k),'.dat'));
if(k==3) this_handle = figure; set(this_handle, 'defaultaxesfontsize', 20); set(this_handle, 'DefaultLineLineWidth', 2); end;

if(k==1) plot(ppot(:,1)/10,ppot(:,2),'Color',[0 0 0],'linewidth',2); end;
if(k==3) plot(ppot(:,1)/10,ppot(:,2),'Color',[0 0 0],'linewidth',2); end;
if(k==4) plot(ppot(:,1)/10,ppot(:,2),'k--','linewidth',2); end;
hold on;
psel=load (strcat(dir,'\psel',num2str(k),'.dat'));
mpt=sum(psel(:,1)==min(psel(:,1)));
npt=length(psel(:,1))/mpt
hold on;
b=1:mpt;
lb=size(b);
A =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0;0.5,0.5,0.5;0.5,0.5,0;1,0.5,0];
A=A(ch,:);
set(gca,'ColorOrder',A);
for j=1:lb(2)
   x(:,j)=psel((b(j)-1)*npt+1:b(j)*npt,1)/10; E(:,j)=psel((b(j)-1)*npt+1:b(j)*npt,2);
end;
B =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0;0.5,0.5,0.5;0.5,0.5,0;1,0.5,0];
A=B(ch,:);
set(gca,'ColorOrder',A);
wa=fa*(E-ones(length(E(:,1)),1)*min(E))+ones(length(E(:,1)),1)*min(E);
   if(k==1) plot(x,wa,'linewidth',2); end;
   if(k==3) plot(x,wa,'linewidth',2); end;
   if(k==4) plot(x,wa,'--','linewidth',2); end;
   hold on;
xlabel('z/nm'); %'Position x [nm]'); %
ylabel('E/eV');  %'Energy [eV]'); %
%for j=1:la(2)
%    plot(psel2((a(j)-1)*203+1:a(j)*203,1),psel2((a(j)-1)*203+1:a(j)*203,2),'--');%'Color',[0.6 0 0]);
%    hold on;
%end;

%this_handle2 = figure;
%    set(this_handle2, 'defaultaxesfontsize', 20);
%    set(this_handle2, 'DefaultLineLineWidth', 2);
%plot(ppot1(:,1),ppot1(:,2));
%axis([-606 0 0 0.8]);
%hold on;
%for j=1:lb(2)
%    plot(psel1((b(j)-1)*203+1:b(j)*203,1),psel1((b(j)-1)*203+1:b(j)*203,2),'Color',[0 0.6 0]);
%    hold on;
%end;