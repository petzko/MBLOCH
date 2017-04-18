classdef TL_solver2 < handle
    %TL_SOLVER Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        dx,dt,Zin,c,modA,modF,V_RF;
        rlgc,bias;
        width,height;
        v_TL,i_TL,i_TLold;
        Bcoeff,Ccoeff,Dcoeff,Ecoeff,Fcoeff;
        N_pts;
        IDX;
        
    end
    
    methods
        function obj = TL_solver2(params,rlgc)
            
        %Definition of constants
        obj.dx = params.dx;              %unit: mm
        obj.dt = params.dt;              %unit: picosecond
        obj.c = params.c;                %unit: mm/picosecond
        obj.IDX = params.IDX;
        obj.N_pts = params.N_pts;
               
        obj.rlgc = rlgc;
        obj.width = params.width;               %unit: mm
        obj.height = params.height;             %unit: mm
        
        obj.bias = params.bias/(obj.height*1e+3); %unit: kV/mm
        obj.modA = params.modA/(obj.height*1e+3); %modulation amplitude unit: kV/mm
        obj.modF = params.modF;          %modulation frequency
        obj.Zin = params.Zin;            %unit: ohm
        obj.rlgc.R = 45*sqrt(obj.modF);  %unit: Ohm/mm --from paper W. Maineult: Microwave modulation of terahertz quantum cascade lasers: a transmission-line approach

        obj.Bcoeff = -obj.width/(obj.height*1e+3)/obj.c/obj.rlgc.C;   % for calculation of v_TL
        
        obj.Ccoeff = 1/(obj.rlgc.C*obj.dx/2/obj.dt+1/2/obj.Zin)/(obj.height*1e+3);
        obj.Dcoeff = (obj.rlgc.C*obj.dx/2/obj.dt-1/2/obj.Zin)*(obj.height*1e+3);
        
        %with consideration of resistance
        obj.Ecoeff = (obj.rlgc.L/obj.dt-1/2*obj.rlgc.R)/(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R); % *I of previous step
        obj.Fcoeff = -1/(obj.dx*(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R))*(obj.height*1e+3/obj.width);% *(V(2)-V(1))
        
        %Initial conditions
        obj.v_TL = obj.bias*ones(obj.N_pts,1);    % V(1)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);            % i(3/2).... i(M+1/2)
        obj.i_TLold = 0; % previous current at node 0
        end
        

        
       function propagate2(obj,J_TL1,J_TL0,t)
        obj.V_RF = obj.modA*sin(2*pi*obj.modF*t)/(obj.height*1e+3);            
%         obj.v_TL(1) = obj.bias+obj.Ccoeff*(obj.Dcoeff*(obj.v_TL(1)-obj.bias)-obj.i_TL(1)*obj.width-obj.width*obj.dx*(J_TL1(1)-J_TL0(1))/2+obj.V_RF/obj.Zin);
        obj.v_TL(1) = (obj.bias+obj.V_RF-(1/2-obj.Zin*obj.rlgc.C*obj.dx/obj.dt)*obj.v_TL(1)-(obj.i_TL(1)-obj.i_TLold)*obj.width*obj.Zin/(obj.height*1e3))/(1/2+obj.Zin*obj.rlgc.C*obj.dx/obj.dt);
        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(2:end-1)-obj.i_TL(1:end-2)+obj.dx*((J_TL1(2:end-1)+J_TL0(2:end-1))/2));
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end)-obj.i_TL(end-1)+(J_TL1(end)+J_TL0(end))*obj.dx/4)*2;

        obj.i_TLold = obj.i_TL(1); 
        obj.i_TL(1:end-1) = obj.i_TL(1:end-1)+obj.Fcoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));
%         obj.i_TL(end) = -obj.i_TL(end-1);


        end

        
          
    end
    
end

