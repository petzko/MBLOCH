close all
clear;

% LOAD a set of simulations 
foldername = 'SIMRES'
fold_cont = dir(foldername)

% load and analyze sim results the name of which follows this regular expression pattern 
pattern = 'test1-mod[^-]\w*' % for T1g = 30ps
% pattern = 'test2-mod[^-]\w*' % for T1g = 20ps
% pattern = 'test3-mod[^-]\w*' % for T1g = 10ps


% cell with all the results 
results = {}

% values of the p parameter
pvals = []
j = 1
for f_idx = 1:length(fold_cont)
   f = fold_cont(f_idx);
   fname =  f.name;
   startidx = regexp(fname,pattern,'match');
   if length(startidx) >0 
       fullname = [foldername,'/',fname]
       % perform the analysis and store the results. 
       results{j} = analyzepulsefunc( fullname )
       pvals(j) = results{j}.p;
       j = j+1;
   end
end

% pretty display of the results
[psorted,pidx] = sort(pvals);
FWHMs = zeros(length(results),1);
FWHMests_alpha0 = zeros(length(results),1); 
FWHMests = zeros(length(results),1); 

T_GR = zeros(length(results),1); 
T1s = zeros(length(results),1);
T2s = zeros(length(results),1);

I0s = zeros(length(results));

POWERs = zeros(length(results),1); 
POWERs_est1 = zeros(length(results),1); 
POWERs_est2 = zeros(length(results),1); 
POWERs_est3 = zeros(length(results),1); 
alphas_est = zeros(length(results),1); 

higher_harmonic_idx = [];
higher_harmonic_ctr = 1; 
for j = 1:length(results)
    if results{pidx(j)}.pulsed == false
        FWHMs(j) = NaN;
        FWHMests_alpha0(j) = NaN;
        FWHMests = NaN;
        
        POWERs_est1(j) = NaN;
        POWERs_est2(j) = NaN;
        POWERs_est3(j) = NaN;
        alphas_est(j) = NaN;
        T_GR(j) =NaN;
    else
        FWHMs(j) = results{pidx(j)}.FWHM; % for a sech pulse
        FWHMests_alpha0(j) = results{pidx(j)}.tau_est*1.763; % for a sech pulse
        FWHMests(j) = results{pidx(j)}.estimates(3)*1.763; % for a sech pulse
        
        POWERs_est1(j) = results{pidx(j)}.P0_est_mW1;
        POWERs_est2(j) = results{pidx(j)}.P0_est_mW2;
        POWERs_est3(j) = results{pidx(j)}.P0_est_mW3;
        
        alphas_est(j) = results{pidx(j)}.estimates(1);

        T_GR(j) = max(results{pidx(j)}.deltaT);
        if results{pidx(j)}.higher_harmonic == true
           higher_harmonic_idx(higher_harmonic_ctr) = j ;
           higher_harmonic_ctr = higher_harmonic_ctr+1;
        end
    end
    I0s(j) = results{pidx(j)}.I0;
    POWERs(j) = results{pidx(j)}.P0_mW;
    T1s(j) = results{pidx(j)}.T1g;
    T2s(j) = results{pidx(j)}.T2g;
    
end

subplot(2,2,1);
ax1 = plot(psorted,FWHMs,psorted,FWHMests_alpha0,psorted,FWHMests);
xlabel('p-val'); ylabel('Pulse FWHM (ps)');
legend('sim value','est. value (\alpha =0)','est value (\alpha \neq 0)');

subplot(2,2,3);
ax2 = plot(psorted,POWERs,psorted,[POWERs_est1,POWERs_est2,POWERs_est3]);
xlabel('p-val'); ylabel('Peak power (mW)');
legend('sim value','est. value 1 (\alpha = 0)','est. value 2 (\alpha = 0)', ... 
    'est. value 3 (\alpha \neq 0)');

ax1(1).Color = [1,0,0];
ax1(2).Color = [1,0,.4];
ax1(2).LineStyle = '--';

ax2(1).Color = [0,0,1];
ax2(2).Color = [.4,0,1];
ax2(2).LineStyle = '--';
ax2(3).LineStyle = '--';
ax2(4).LineStyle = '--';

% in case the laser mode locks in higher harmonic regime (2 pulses per
% roudn trip etc. as in the case of T1g = 10 ps )
% place a cross sign everywhere on the plots!  

for idx = 1:length(higher_harmonic_idx)
    h_idx = higher_harmonic_idx(idx);
    subplot(2,2,1);
    hold on
    plot(psorted(h_idx),FWHMs(h_idx),'bx');
    plot(psorted(h_idx),FWHMests_alpha0(h_idx),'bx');
    subplot(2,2,3);
    hold on
    plot(psorted(h_idx),POWERs(h_idx),'rx');
    plot(psorted(h_idx),POWERs_est1(h_idx),'rx'); 
     plot(psorted(h_idx),alphas_est(h_idx),'rx'); 
end
subplot(2,2,[2,4]);
plot(psorted,T_GR,'-k'); 
xlabel('p-val');
ylabel('Gain recovery time (ps)');



I0guess = (psorted.'-1).^2./(T1s.*T2s);


figsize = [0,0,0.4,0.2]
figure('units','normalized','position',figsize)
subplot(1,2,1);
ax = plotyy(psorted,FWHMs,psorted,POWERs);
set(ax(1).YLabel,'String','Pulse FWHM (ps)');
ax(1).XLim=[psorted(1),psorted(end)];
ax(1).Children.Color = [1,0,0];
ax(1).YColor = [1,0,0];

ax(2).XLim=[psorted(1),psorted(end)];
ax(2).Children.Color = [0,0,1];
ax(2).YColor = [0,0,1];



set(ax(1).XLabel,'String','p-val');
set(ax(2).YLabel,'String','Peak power (mW)');

fh = subplot(1,2,2);
ax = plot(psorted,T_GR,'-k'); 
xlim([psorted(1),psorted(end)]);
fh.YAxisLocation = 'right';

xlabel('p-val');
ylabel('Gain recovery time (ps)');

