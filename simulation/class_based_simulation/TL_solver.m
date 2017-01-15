classdef TL_solver < handle
    %TL_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dx,dt,Zin,Vs,c;
        rlgc,current,bias;
        width,height;
        v_TL,i_TL,J_TL;
        Acoeff, Bcoeff,Ccoeff,Dcoeff;
        N_pts;
        IDX;
        
    end
    
    methods
        function obj = TL_solver(params,rlgc)
            
        %Definition of constants
        obj.dx = params.dx;              %unit: mm
        obj.dt = params.dt;              %unit: picosecond
        obj.c = params.c;                %unit: mm/picosecond
        obj.IDX = params.IDX;
        obj.N_pts = params.N_pts;
        
        obj.Zin = params.Zin*1e-3;       %unit: kV/A
        obj.bias = params.bias/10;       %unit: kV/mm
        obj.Vs = params.Vs*1e-3;         %unit: kV -->initial voltage
        obj.current = params.current;    %unit:A/mm -->initial current
        
        obj.rlgc = rlgc;
        obj.width = params.width;               %unit: mm
        obj.height = params.height;             %unit: mm

        obj.Acoeff = -obj.height/obj.width/obj.c/obj.rlgc.L*1e+3;   % for calculation of i_TL  % ???
        obj.Bcoeff = -obj.width/obj.height/obj.c/obj.rlgc.C*1e-3;   % for calculation of v_TL  % ???
%         obj.Ccoeff = obj.rlgc.C*obj.dx/2/obj.dt+1/2/(obj.Zin*1e+3);
%         obj.Dcoeff = 1/2/(obj.Zin*1e+3)-obj.rlgc.C*obj.dx/2/obj.dt;
        
        %Initial conditions
        obj.v_TL = obj.bias*ones(obj.N_pts,1);    % V(1)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);          % I(0)  I(1/2)  I(3/2).... I(M+1/2)

        obj.i_TL(1) = obj.current;
        for mm = 2: obj.N_pts-1
        obj.i_TL(mm) = obj.i_TL(1)*(obj.N_pts-mm)/obj.N_pts;
        end

        end
        
        function propagate_1(obj,V2_TL,J_TL)
            
        obj.J_TL = J_TL;

        obj.v_TL(1) = obj.Vs/obj.height-V2_TL;
%         obj.v_TL(1) = obj.v_TL(1)+2*obj.Bcoeff*(obj.i_TL(2)-obj.i_TL(1)+obj.dx*obj.J_TL(1)/2);
        obj.v_TL(2:end) = obj.v_TL(2:end)+obj.Bcoeff*(obj.i_TL(2:end)-obj.i_TL(1:end-1)+obj.dx*obj.J_TL(2:end));

%         obj.i_TL(1) = obj.current; obj.i_TL(end) = -obj.i_TL(end-1);
        obj.i_TL(1:end-1) = obj.i_TL(1:end-1)+obj.Acoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));


        end
       
        function propagate_2(obj,I1_TL,J_TL)
            
        obj.J_TL = J_TL;

        obj.v_TL(1) = obj.v_TL(1)+2*obj.Bcoeff*(obj.i_TL(1)-I1_TL+obj.dx*obj.J_TL(1)/2);
        
        obj.v_TL(2:end) = obj.v_TL(2:end)+obj.Bcoeff*(obj.i_TL(2:end)-obj.i_TL(1:end-1)+obj.dx*obj.J_TL(2:end));
        obj.i_TL(1:end-1) = obj.i_TL(1:end-1)+obj.Acoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));


        end

        
          
    end
    
end

