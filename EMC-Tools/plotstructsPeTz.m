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
	plotstructsPeTz(dir,ch)
%e.g., ch=[3 4 1 5 2] means that the level with the second-lowest energy
%should be colored with color '4', which might indicate the lower laser
%level, and so on

fa=.5;
ppot=load (strcat(dir,'\ppot1.dat'));
plot(ppot(:,1),ppot(:,2),'Color',[0 0 0],'linewidth',2);

hold on;
% psel=load (strcat(dir,'\psel',num2str(k),'.dat'));
WF_grid =  load(strcat(dir,'\psi.dat'));
Ens = load('E_BOUND');  
NrWfs = length(Ens(:,1))/4;
Npts = length(WF_grid(:,1))/(4*NrWfs); 


B =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0;0.5,0.5,0.5;0.5,0.5,0;1,0.5,0];
A=B(ch,:);
set(gca,'ColorOrder',A);

for i = 1:4*NrWfs
    Phi = WF_grid((i-1)*Npts+1:i*Npts,2);
    x = WF_grid((i-1)*Npts+1:i*Npts,1);
    plot(x,fa*abs(Phi).^2+Ens(i,2),'Linewidth',2.0);
end




xlabel('z/nm'); %'Position x [nm]'); %
ylabel('E/eV');  %'Energy [eV]'); %

