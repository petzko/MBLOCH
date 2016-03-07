

fold = dir(pwd);
phases  = {} ;
Ymods = {};
el_ctr =1; 
for fil =1:length(fold)
    if(findstr(fold(fil).name,'Phase_const'))
        phases{el_ctr} = fold(fil).name;
        el_ctr = el_ctr+ 1;
    end
end
c =   0.0833; dx= 0.0013;
N2 = length(Ymod);
k_ = (2*pi/N2)*[-N2/2:N2/2-1].'*1/dx;
dfigure; hold on; 
for phi = 1:length(phases)
    load(phases{phi});
      
    plot((k_*c/2/pi),fftshift(angle(Ymod)));
end
legend(phases)