function testWFplots(direc)

files = dir(direc);
for i = 1:length(files)
    if(findstr(files(i).name,'ENs_at'))
        ENfile = [direc '\' files(i).name]; 
    else if (findstr(files(i).name,'POT_at'))
            POTfile =[direc '\' files(i).name]; 
        else if (findstr(files(i).name,'WFs_at'))
            WFfile = [direc '\' files(i).name ]; 
            end
        end
    end

end


Ens = load(ENfile);  
NrWfs = length(Ens(:,1))/4; 
fopen(WFfile,'r');

WF_grid = dlmread(WFfile,'\t',50,0);
POT_grid = load(POTfile); 

hold on ; 
plot(POT_grid(:,1),POT_grid(:,2))
Npts = length(WF_grid(:,1))/(4*NrWfs); 

for i = 1:4*NrWfs
    
    Phi = WF_grid((i-1)*Npts+1:i*Npts,2);
    x = WF_grid((i-1)*Npts+1:i*Npts,1);
    plot(x,abs(Phi).^2+Ens(i,2));
    

    
end




