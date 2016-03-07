%%%% extract the IV characteristic of the laser
clear all; 
clc; 
% TEMPLVLS = 10:10:200;
TEMPLVLS = 50;

ntemps = length(TEMPLVLS);

W_sf = zeros(ntemps,4*4); 
W_sf_op = W_sf; W_sf_ac = W_sf; W_sf_ee = W_sf; W_sf_imp = W_sf; W_sf_if = W_sf; 


Energies  = zeros(ntemps,4); 

maindir = pwd; 

%current density in A/m^2
    global sf; global sf_op;  global sf_ac; global sf_ee; global sf_imp; global sf_if;     
    
    % extract  the energies in the following order Injector -> ULL -> LLL -> DEPOP
    lvls = [13 12 10 8];
    for vidx = 1:ntemps  
            folder_ = [maindir];
            %get the current energies.. 
            Ens = load([folder_ '\E_BOUND']);
            [sortedEns, order]  = sort(Ens(1:end,2));
%             Energy_lvls = sortedEns(lvls,1);
            Energies(vidx,:) = sortedEns(lvls);
    end
     
for i=1:length(lvls)
     for j = 1:length(lvls) 
        for vidx = 1:ntemps  
            folder_ = [maindir];
            scan(folder_,lvls(i), lvls(j), 1:5,[],3);
            W_sf(vidx,j+length(lvls)*(i-1)) = sf; W_sf_op(vidx,j+length(lvls)*(i-1)) = sf_op; W_sf_ac(vidx,j+length(lvls)*(i-1)) = sf_ac; 
            W_sf_ee(vidx,j+length(lvls)*(i-1)) = sf_ee; W_sf_imp(vidx,j+length(lvls)*(i-1)) = sf_imp; W_sf_if(vidx,j+length(lvls)*(i-1))= sf_if; 
        end
     end
end

W = zeros(length(lvls),length(lvls));
INJ = 1; ULL = 2; LL = 3; DEPOP = 4; 
IDX = [INJ ULL LL DEPOP];
for i=IDX
     for j=IDX
        W(i,j) =  W_sf(1,j+length(lvls)*(i-1)); 
     end
end
W =W/1E12



% filename = 'temp_rates.dat';
% fileid = fopen(filename,'w');
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% 
% fprintf(fileid,'	Energies (Inj , ULL, LLL, DEPOP) -> \r\n');
% fprintf(fileid,' TEMP |\r\n');
% fprintf(fileid,'      ---------------------------------------------------------------\r\n');
% 
% 
% 
% for vidx = 1:ntemps
%     fprintf(fileid,'%3i K    ; \t %.5f ; \t %.5f ; \t %.5f ; \t %.5f ; \r\n',vidx*10,Energies(vidx,1),Energies(vidx,2),Energies(vidx,3),Energies(vidx,4));
% end
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% 
% fprintf(fileid,' Total scattering rates (THz): \r\n');
% fprintf(fileid,' TEMP ; Transition (Inj -> ULL, Inj -> LLL, ...., DEPOP-> DEPOP)\r\n');
% 
% for tmpidx = 1:ntemps
%     fprintf(fileid, '%3iK ; ',tmpidx*10);
%     for i = 1:length(lvls) 
%         for j = 1:length(lvls) 
%              fprintf(fileid,' %.5f ;',W_sf(tmpidx,j+length(lvls)*(i-1))/1E12); 
%         end
%     end
%    fprintf(fileid,'\r\n');        
% end
% 
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% 
% fprintf(fileid,' Optical phonon scattering rates (THz):\r\n');
% fprintf(fileid,' TEMP ; Transition (Inj -> ULL, Inj -> LLL, ...., DEPOP-> DEPOP)\r\n');
% 
% for tmpidx = 1:ntemps
%     fprintf(fileid, '%3iK ; ',tmpidx*10);
%     for i = 1:length(lvls) 
%         for j = 1:length(lvls) 
%              fprintf(fileid,' %.5f ; ',W_sf_op(tmpidx,j+length(lvls)*(i-1))/1E12); 
%         end
%     end
%    fprintf(fileid,'\r\n');        
% end
% 
% 
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% 
% fprintf(fileid,' Acoustic phonon scattering rates (THz):\r\n');
% fprintf(fileid,' TEMP ; Transition (Inj -> ULL, Inj -> LLL, ...., DEPOP-> DEPOP)\r\n');
% 
% for tmpidx = 1:ntemps
%     fprintf(fileid, '%3iK ; ',tmpidx*10);
%     for i = 1:length(lvls) 
%         for j = 1:length(lvls) 
%              fprintf(fileid,' %.5f ;',W_sf_ac(tmpidx,j+length(lvls)*(i-1))/1E12); 
%         end
%     end
%    fprintf(fileid,'\r\n');        
% end
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% fprintf(fileid,' Electron-electron scattering rates (THz) : \r\n');
% fprintf(fileid,' TEMP ; Transition (Inj -> ULL, Inj -> LLL, ...., DEPOP-> DEPOP)\r\n');
% 
% for tmpidx = 1:ntemps
%     fprintf(fileid, '%3iK ; ',tmpidx*10);
%     for i = 1:length(lvls) 
%         for j = 1:length(lvls) 
%              fprintf(fileid,' %.5f ;',W_sf_ee(tmpidx,j+length(lvls)*(i-1))/1E12); 
%         end
%     end
%    fprintf(fileid,'\r\n');        
% end
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% fprintf(fileid,' Interface-roughness scattering rates (THz): \r\n');
% fprintf(fileid,' TEMP ; Transition (Inj -> ULL, Inj -> LLL, ...., DEPOP-> DEPOP)\r\n');
% 
% for tmpidx = 1:ntemps
%     fprintf(fileid, '%3iK ; ',tmpidx*10);
%     for i = 1:length(lvls) 
%         for j = 1:length(lvls) 
%              fprintf(fileid,' %.5f ;',W_sf_if(tmpidx,j+length(lvls)*(i-1))/1E12); 
%         end
%     end
%    fprintf(fileid,'\r\n');        
% end
% 
% 
% fprintf(fileid,'\r\n##########################################################################\r\n');
% fprintf(fileid,' Impurity scattering rates (THz): \r\n');
% fprintf(fileid,' TEMP ; Transition (Inj -> ULL, Inj -> LLL, ...., DEPOP-> DEPOP)\r\n');
% 
% for tmpidx = 1:ntemps
%     fprintf(fileid, '%3iK ; ',tmpidx*10);
%     for i = 1:length(lvls) 
%         for j = 1:length(lvls) 
%              fprintf(fileid,' %.5f ;',W_sf_imp(tmpidx,j+length(lvls)*(i-1))/1E12); 
%         end
%     end
%    fprintf(fileid,'\r\n');        
% end
% 
% fclose(fileid); 
