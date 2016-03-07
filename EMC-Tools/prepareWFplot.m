function [ output_args ] = prepareWFplot(direc)
%PREPAREWFPLOT Summary of this function goes here
%   Detailed explanation goes here

Ens = load([direc '\E_BOUND']); 
Ens = Ens(:,2);
Wfs = dlmread([direc '\WFs_at_10.2_kVpcm(07.08.2015).txt'],['\t'],50, 0);

x = Wfs(1:end,1); WF = Wfs(1:end,2); 


NWF = length(Ens);

Nx = length(WF)/NWF;
xinterp = linspace(x(1),x(Nx),abs(x(1))).';
x_all = []; WF_all = [];
x1 = x(1:Nx);
hold on;
ppot=load (strcat(direc,'\ppot1.dat'));
plot(ppot(:,1),ppot(:,2),'Color',[0 0 0],'linewidth',2); 

for i = 1:NWF
    x_all =[x_all ; xinterp];
    plot(x1,abs(WF((i-1)*Nx+1:i*Nx)).^2+Ens(i))
    WF((i-1)*Nx+1:i*Nx) =  abs(WF((i-1)*Nx+1:i*Nx)).^2+Ens(i);
    wf = interp1(x((i-1)*Nx+1:i*Nx),WF((i-1)*Nx+1:i*Nx),xinterp);
    WF_all = [WF_all; wf];
end

A = [x_all  WF_all];


fileId = fopen('wfs2plot.dat','w');
fprintf(fileId, '       %i  %.15E \r\n',A');
fclose(fileId);

end

