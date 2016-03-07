clear;
% close all; 
% figure;
% folder = '11p2-disp-strong';
folder = '11-disp-standard'
fileTA = dir([ pwd '/' folder] ); 


allctr =1; 
amplis = [];
b1s = [] ;
b2s = [] ;
b3s = [];
nam = '.mat';
amplitudes_vec = []; 
beta_1_vec = [];
beta_2_vec = []; close 
beta_3_vec = []; 

hfig = dfigure('DName','Dispersion analysis'); 
subplot(1,2,1);  hold on; 
for i = 1:length(fileTA)
    fname =[pwd '/' folder '/' fileTA(i).name]; 
    if(strfind(fname,nam))
        load(fname); 
        dispanalysis_IV;
        plot(f_+f0,n_Re);
        legInf{allctr} = [num2str(ampl)];
        amplitudes_vec(allctr) = ampl; 
        beta_1_vec(allctr) = b1; 
        beta_2_vec(allctr) = b2; 
        beta_3_vec(allctr) = b3; 
        
        allctr = allctr + 1; 
    end
end
% dlegend(legInf,'Initial Amplitude:');
dlegend(legInf,'Init. amplitude (2\pi\timesTHz):');
xlabel('Freq. (THz)');
ylabel('(n_{Re})'); xlim([3.4,4.4])


subplot(2,2,2); hold on;
allctr =1; 
for i = 1:length(fileTA)
    fname =[pwd '/' folder '/' fileTA(i).name]; 
    if(strfind(fname,nam))
        load(fname); 
        plotON = false;
        dispanalysis_IV;
        plot(f_+f0,gn*10);
        legInf{allctr} = [num2str(ampl)]; 
        allctr = allctr + 1; 
    end
end
xlabel('Freq. (THz)');xlim([3.4,4.4])
ylabel('Gain (1/cm)');
set(hfig,'position', [100, 100, 500, 200]) 


subplot(2,2,[3 4]);
for i = 1:length(fileTA)
    fname =[pwd '/' folder '/' fileTA(i).name]; 
    if(strfind(fname,nam))
        load(fname); 
        plotON = false;
        dispanalysis_IV;
        plot(f_+f0,gn*10);
        legInf{allctr} = [num2str(ampl)]; 
        allctr = allctr + 1; 
    end
end
