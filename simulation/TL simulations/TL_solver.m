classdef TL_solver < handle
    %TL_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dx,dt,Zin,Vs,c,modA,modF,Vs2,Res_TL;
        rlgc,current,bias,k;
        width,height;
        v_TL,i_TL;
        Acoeff, Bcoeff,Ccoeff,Dcoeff,Ecoeff,Fcoeff;
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
               
        obj.rlgc = rlgc;
        obj.width = params.width;               %unit: mm
        obj.height = params.height;             %unit: mm
        
        obj.bias = params.bias/(obj.height*1e+3); %unit: kV/mm
        obj.Vs = params.Vs*1e-3;         %unit: kV -->initial voltage
        obj.modA = params.modA/(obj.height*1e+3); %modulation amplitude unit: kV/mm
        obj.modF = params.modF;          %modulation frequency
        obj.Zin = params.Zin;            %unit: ohm
        obj.Res_TL = 45*sqrt(obj.modF)*params.length; % QCL resistance unit: Ohm --from paper W. Maineult: Microwave modulation of terahertz quantum cascade lasers: a transmission-line approach
        obj.k = obj.Res_TL/(obj.Zin+obj.Res_TL);  % scale factor, ratio of AC current and DC 
%         obj.k = 0.04;

        obj.Acoeff = obj.k*obj.Zin*obj.rlgc.C*obj.dx+1;   % New boundary condition with regard to wire impedance
        obj.Bcoeff = -obj.width/(obj.height*1e+3)/obj.c/obj.rlgc.C;   % for calculation of v_TL
        
        obj.Ccoeff = 1/(obj.rlgc.C*obj.dx/2/obj.dt+1/2/obj.Zin)/(obj.height*1e+3);
        obj.Dcoeff = (obj.rlgc.C*obj.dx/2/obj.dt-1/2/obj.Zin)*(obj.height*1e+3);
        
        %with consideration of resistance
        obj.Ecoeff = (obj.rlgc.L/obj.dt-1/2*obj.rlgc.R*obj.k)/(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R*obj.k); % *I(old)
        obj.Fcoeff = -1/(obj.dx*(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R*obj.k))*(obj.height*1e+3/obj.width);% *(V(2)-V(1))
        
        %Initial conditions
        obj.v_TL = obj.bias*ones(obj.N_pts,1);    % V(1)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);            % I(3/2).... I(M+1/2)

        end
        

        
       function propagate(obj,J_TL1,t)
        obj.Vs2 = obj.bias+obj.modA*obj.k*sin(2*pi*obj.modF*t);            
%         obj.v_TL(1) = obj.Ccoeff*(obj.Dcoeff*obj.v_TL(1)-obj.i_TL(1)*obj.width-obj.width*obj.dx*J_TL1(1)/2+obj.Vs2/obj.Zin);
%         obj.v_TL(1) = (obj.Acoeff-2)/obj.Acoeff*obj.v_TL(1)+2/obj.Acoeff*(obj.Vs2-obj.k*obj.Zin*obj.width*(obj.dx/2*J_TL1(1)+obj.i_TL(1))*1e-3/obj.height);
%        obj.Vs2 = obj.bias*(1+obj.modA*obj.k*sin(2*pi*obj.modF*t));
        obj.v_TL(1) = obj.Vs2;
        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(2:end-1)-obj.i_TL(1:end-2)+obj.dx*J_TL1(2:end-1));
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end)-obj.i_TL(end-1)+J_TL1(end)*obj.dx/2)*2;

        obj.i_TL(1:end-1) = obj.Ecoeff*obj.i_TL(1:end-1)+obj.Fcoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));
        obj.i_TL(end) = 0;


        end

        
          
    end
    
end

