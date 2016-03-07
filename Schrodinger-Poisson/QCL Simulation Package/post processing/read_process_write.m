function read_process_write(folder)

files = dir(folder);

for i = 1:length(files)-1
    clear('Energies','xWF','wFunctions','x','V');
    fname = files(i).name;
    if(strfind(fname,'WF'))  
    
        fid = fopen([folder '\' fname],'r'); 
        %get first line -> number of WFs 
        tline = fgetl(fid);
        NrWF = str2num(tline); 
        
        %get second line -> number of grid points
        tline = fgetl(fid);
        NrPts = str2num(tline); 
  
        %get third line -> Biasbase
        tline = fgetl(fid);
        Biasbase = str2num(tline); 
        
        %skip line
        tline = fgetl(fid);
        %get fifth line -> energies! 
        tline = fgetl(fid);
        Energies = str2num(tline);
        
        %skip line
        tline = fgetl(fid);
        tline = fgetl(fid);

        ctr = 1; 
        A = zeros(NrPts,NrWF+1) ;

        while (ischar(tline))
           A(ctr,:) = str2num(tline) ; 
           tline = fgetl(fid);
           ctr  = ctr+1;
        end

        fclose(fid);

        xWF = A(:,1)/10;
        wFunctions = A(:,2:end);
      
        TB2EB_transform; 
        MC_PRINT
    end
end