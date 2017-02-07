classdef TL_solver < handle
    %TL_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dx,dt,Zin,Vs,c;
        rlgc,current,bias;
        width,height;
        v_TL,i_TL;
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
        
        obj.bias = params.bias/10;       %unit: kV/mm
        obj.Vs = params.Vs*1e-3;         %unit: kV -->initial voltage
        obj.current = params.current;    %unit:A/mm -->initial current
        obj.Zin = params.Zin;            %unit: ohm
        
        obj.rlgc = rlgc;
        obj.width = params.width;               %unit: mm
        obj.height = params.height;             %unit: mm

        obj.Acoeff = -(obj.height*1e+3)/obj.width/obj.c/obj.rlgc.L;   % for calculation of i_TL
        obj.Bcoeff = -obj.width/(obj.height*1e+3)/obj.c/obj.rlgc.C;   % for calculation of v_TL
        
        obj.Ccoeff = 1/(obj.rlgc.C*obj.dx/2/obj.dt+1/2/obj.Zin)/(obj.height*1e+3);
        obj.Dcoeff = (obj.rlgc.C*obj.dx/2/obj.dt-1/2/obj.Zin)*(obj.height*1e+3);
        
        %Initial conditions
        obj.v_TL = obj.bias*ones(obj.N_pts,1);    % V(1)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);            % I(3/2).... I(M+1/2)

        for mm = 1: obj.N_pts-1
        obj.i_TL(mm) = obj.current*(obj.N_pts-mm)/obj.N_pts;
        end
        obj.i_TL(end)= -obj.i_TL(end-1);

        end
        
        function propagate_1(obj,J_TL1)
            
% %         i1new = obj.i_TL(1)+J_TL1(1)*obj.dx;
%         obj.v_TL(1) = obj.Vs/obj.height;
        obj.v_TL(1) = obj.Ccoeff*(obj.Dcoeff*obj.v_TL(1)-obj.i_TL(1)*obj.width-obj.width*obj.dx*J_TL1(1)/2+obj.Vs*1e+3/obj.Zin);
% %         obj.v_TL(1) = (obj.Vs - (i1new-i1old)*obj.width*(obj.Zin*1e-3))/obj.height;
        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(2:end-1)-obj.i_TL(1:end-2)+obj.dx*J_TL1(2:end-1));
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end)-obj.i_TL(end-1)+J_TL1(end)*obj.dx/2)*2;
%         obj.i_TL(1) = obj.current; obj.i_TL(end) = -obj.i_TL(end-1);
        obj.i_TL(1:end-1) = obj.i_TL(1:end-1)+obj.Acoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));
        obj.i_TL(end) = -obj.i_TL(end-1);


        end
       
        function propagate_2(obj,J1_tot,J_TL2)
            

        obj.v_TL(1) = obj.v_TL(1)+obj.Bcoeff*(obj.i_TL(1)-(J1_tot-obj.v_TL(1)*1e3*obj.height/obj.Zin/obj.width)+J_TL2(1)*obj.dx/2)*2;
% %         obj.v_TL(1) = obj.Ccoeff*(obj.Dcoeff*obj.v_TL(1)-obj.i_TL(1)*obj.width-obj.width*obj.dx*J_TL2(1)/2+J1_tot*obj.width);
% % %         obj.v_TL(1) = obj.v_TL(1)+obj.Bcoeff*(obj.i_TL(1)-(J1_tot-obj.v_TL(1)*obj.height*1e3/obj.Zin/obj.width)+obj.dx*J_TL2(1)/2);

        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(2:end-1)-obj.i_TL(1:end-2)+obj.dx*J_TL2(2:end-1));
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end)-obj.i_TL(end-1)+J_TL2(end)*obj.dx/2)*2;
%         obj.i_TL(1) = I1_TL-J_TL2(1)*obj.dx/2;
        obj.i_TL(1:end-1) = obj.i_TL(1:end-1)+obj.Acoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));
        obj.i_TL(end) = -obj.i_TL(end-1);        


        end

        
          
    end
    
end

