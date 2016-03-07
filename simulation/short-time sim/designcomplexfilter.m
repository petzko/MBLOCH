% n = 50;     % Group delay of a linear phase filter would be 25.
% gd = 12;    % Set the desired group delay for the filter.
% f1=linspace(0,.25,30); % Define the first stopband frequencies.
% f2=linspace(.3,.56,40);% Define the passband frequencies.
% f3=linspace(.62,1,30); % Define the second stopband frequencies.
% h1 = zeros(size(f1));  % Specify the filter response at the freqs in f1.
% h2 = exp(-1j*pi*gd*f2); % Specify the filter response at the freqs in f2.
% h3 = zeros(size(f3));  % Specify the response at the freqs in f3.
% d=fdesign.arbmagnphase('n,b,f,h',50,3,f1,h1,f2,h2,f3,h3);
% D = design(d,'equiripple');
% fvtool(D,'Analysis','freq');
% 
% gd = 12;    % Set the desired group delay for the filter.
clear;
Nf = 1000; 

F1=linspace(-1,1,Nf).'; % Define the first stopband frequencies.
H1 = hanning(Nf).*exp(1i*F1.^2);  % Specify the filter response at the freqs in f1.
f = fdesign.arbmagnphase('N,F,H',100,F1,H1);
Hd = design(f,'equiripple');
hfvt = fvtool(Hd, 'Color','w');

hfvt(2) = fvtool(Hd,'Analysis','phase','Color','white');

ax = hfvt(2).CurrentAxes;
ax.NextPlot = 'add';
pidx = find(F1>=0);
plot(ax,F1,[fliplr(unwrap(angle(H1(pidx-1:-1:1)))) ... % Mask
    unwrap(angle(H1(pidx:end)))],'k--')