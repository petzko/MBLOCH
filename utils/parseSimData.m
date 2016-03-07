function settings = parseSimData(filename,settings)

assert(exist(filename) == 2,'the specified file does not exist. please try again'); 

optnames = {'name','zUL','Energies','O13','INJ','ULL','LLL','DEPOP','HTB','Wmtx'};
n = length(optnames);
comments = '#%!';
optvals = cell(n,1);

for i = 1:n
    optvals{i} = findoption(filename,optnames{i},comments);
end


for i = 1:n
    settings.(optnames{i}) = optvals{i};
end


    function val = findoption(filename,optname,comments)
        
        fid = fopen(filename,'r'); 
        
        tline = fgetl(fid);
        val = 'Unknown'; 
        
        while ischar(tline)
            skip = false;
            if(length(tline) >= 1)
                for idx = 1:length(comments)
                    %skip comments 
                        if(tline(1) == comments(idx))
                         skip = true; break;
                        end
                end
                
                if (skip) tline = fgetl(fid); continue; end;
                    % parse option
                    
                tmp = strread(tline,'%s');
                %verify for empty strings. 
                if(length(tmp) <=2)
                    tline = fgetl(fid);
                    continue;
                end
                
                name = tmp(1); 
                
                tmpval = tmp(3:end); 
                strval = [];
                for idx=1:length(tmpval)
                    strval = [strval ' ' tmpval{idx}]; 
                end
                
                
                if ~strcmp(name,optname)
                    tline = fgetl(fid);
                    continue; 
                end

                %found our options 
                display(['Match: ' name optname]);
            
                if(strcmp(name,'name'))
                    val = strtrim(strval); %% remove leading and trailing white spaces
                else
                    val = str2num(strval);
                end
                
                break;   
               end
            tline = fgetl(fid);
        end
        fclose(fid); 
        
        
    end

end

