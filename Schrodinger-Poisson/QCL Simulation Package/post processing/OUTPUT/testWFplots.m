function testWFplots(direc,wfreadOFFSET)

files = dir(direc);

%find the wavefunctions, potential profile and EN files. 
flag = -1; 
for i = 1:length(files)
    file = files(i).name;
    if(strfind(file,'ENs'))
            ENfile = [direc '\' file];
            flag = flag*-1;
    else if strfind(file,'WFs')
           WFfile = [direc '\' file]; 
            flag = flag*-1;
        else if strfind(file,'POT')
            POTfile =[direc '\' file]; 
            flag = flag*-1;
        end
        end
    end
end

%could not find files
if(flag < 0 )
    exit('Folder does not contain either on or more of the neccessary files. Aborting');
end


%load the energies 
Ens = load(ENfile);  
NrWfs = length(Ens(:,1))/4; 

%read the wavefunctions! 

WF_grid = dlmread(WFfile,'\t',wfreadOFFSET,0);
POT_grid = load(POTfile); 

hold on ; 
plot(POT_grid(:,1),POT_grid(:,2))
Npts = length(WF_grid(:,1))/(4*NrWfs); 

for i = 1:4*NrWfs
    
    Phi = WF_grid((i-1)*Npts+1:i*Npts,2);
    x = WF_grid((i-1)*Npts+1:i*Npts,1);
    plot(x,0.05*abs(Phi).^2/max(abs(Phi).^2)+Ens(i,2));
    

    
end
