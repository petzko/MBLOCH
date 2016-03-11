mkdir(['OUTPUT/' num2str(Biasbase) ])
PSIname = ['OUTPUT/' num2str(Biasbase) '/WFs_at_' num2str(Biasbase) '_kVpcm(07.08.2015).txt']; 
ENname = ['OUTPUT/' num2str(Biasbase) '/ENs_at_' num2str(Biasbase) '_kVpcm(07.08.2015).txt'];
Potname = ['OUTPUT/' num2str(Biasbase) '/POT_at_' num2str(Biasbase) '_kVpcm(07.08.2015).txt'];

x_all = [];  wf_all = []; E_all = []; 
skipsize = 3;  dz = 0.3;

x_A = [x_fin(1):0.1:x_fin(end)].';
TB_A = zeros(length(x_A),NrWF,Nper);
V_A = interp1(x_fin,V_fin,x_A,'linear',0);

for p = 1:Nper
    for i = 1:NrWF
    TB_A(:,i,p) = interp1(x_fin,TB_WF2(:,i,p),x_A,'linear',0);
    end
end

x_neg = x_A - x_A(end);
dz = (x_neg(2)- x_neg(1))*skipsize;
data = [round(x_neg*10) , V_A];

figure;  plot(x_neg,V_A,'-k','Linewidth',2.0);


hold on

selector = {[3,4,5], [1,2,3,4,5] , [1,2,3,4,5] , [1,2,3,4,5] ,[1,2]};
B =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0;0.5,0.5,0.5;0.5,0.5,0;1,0.5,0];
colors=B(1:NrWF,:);
% set(gca,'ColorOrder',A);

for p = 1:length(selector)
    wfs = selector{p};
    for i = 1:length(wfs)
        E_all = [E_all Ens(wfs(i),p)]; 
        x_all = [x_all ; x_neg(1:skipsize:end) ];
        wf_all = [wf_all ; TB_A(1:skipsize:end,wfs(i),p) ];
%         plot(x_neg(1:skipsize:end),abs(TB_A(1:skipsize:end,wfs(i),p)).^2 + Ens(wfs(i),p),'Linewidth',2.0,'color',colors(wfs(i),:));
    end
end


A = [round(x_all*10)  wf_all];


%%% do the printing!!! 
potfileId = fopen(Potname,'w');
fprintf(potfileId,'\t %d %f \r\n',data');
fclose(potfileId);

fileId = fopen(PSIname,'w'); 
printEnergies = [E_all ;  zeros(1, length(E_all)) ; zeros(1, length(E_all));ones(1, length(E_all))];

someIntegers = [20 6 5  4  3  3   11   12  13   14  15   6  7   8  9    10    11   12   13    14  15   6  7  8  9  10];
fprintf(fileId,'         %4d \r\n',round(length(x_neg(1:skipsize:end))));
fprintf(fileId,'    %.14f \r\n',round(dz*10));
fprintf(fileId,'        %d \r\n', round(x_neg(1)*10));
fprintf(fileId,'   %.15E \r\n',0.0);
fprintf(fileId,'           %2d \r\n',round(someIntegers));

fprintf(fileId,'    %.10f         %.9f         %.9f         %.9f\r\n',printEnergies);
fprintf(fileId, '       %i  %.15E \r\n',A');
fclose(fileId); 

fileId = fopen(ENname,'w'); 
fprintf(fileId,'           %2d  %.10f \n',[ [1:length(E_all)] ; E_all]);
fclose(fileId); 

savefile = ['OUTPUT/' num2str(Biasbase)  '/CHECKPOINT_at' num2str(Biasbase) '_kV_cm'];
save(savefile);
