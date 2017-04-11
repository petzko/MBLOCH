function res = varfunc(varargin)
    if length(varargin)>0
        name = varargin{1}
        val = varargin{2}
        display( [name, ' ', val] ) 
    end
end