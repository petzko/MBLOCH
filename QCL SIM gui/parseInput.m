function settings = parseInput(filename)

assert(exist(filename) == 2,'the specified file does not exist. please try again'); 

optnames = {'out_coupling','shb','tch','lch','N','D','Ltot','n','z23','E1','E2','E3',...
    'O13','INJ','ULL','LLL','DEPOP','W_inj','W_ull','W_lll','W_dep','dN','Overlap','loss','simRT','plotCtr','recordRT'};
n = length(optnames);
comments = '#%!';

for i = 1:n
    optvals{i} = findoption(filename,optnames{i},comments);
end


settings = struct();
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
                val = str2num(strval);
                break;   
               end
            tline = fgetl(fid);
        end
        fclose(fid); 
        
        
    end

end

