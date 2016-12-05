function settings = input_parser(filename)

assert(exist(filename) == 2,'the specified file does not exist. please try again');

optnames = {'name','zUL','E_UL','E_LL','T_1','T_2','L', 'R', ...
            'Gamma','Overlap','nTHz','nRF','tch','lch','rho_u_0','rho_l_0',...
            'linear_loss','Ncarriers_cm','H','W','Grnd_','Ex_','Spin_', ...
            'D','Lp','Ld','dN','modA','modF','bias','voltage','current','Ltot', ...
            'N_pts','shb','linear_loss','Tdeph_13','Tdeph_12','Tdeph_32','Zin',...
            'width','height','length','Vs'};


n = length(optnames);
comments = '#%!';
found_options = {};
found_idx = 1;
found_vals = {};
for i = 1:length(optnames)
    val = findoption(filename,optnames{i},comments);
    if ~strcmp(val,'Unknown')
        found_options{found_idx} = optnames{i};
        found_vals{found_idx} = val;
        found_idx = found_idx + 1;
    end
end


for i = 1:length(found_options)
    settings.(found_options{i}) = found_vals{i};
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

