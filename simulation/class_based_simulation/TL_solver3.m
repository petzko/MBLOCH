classdef TL_solver3 < handle
    %TL_SOLVER Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        dx,dt,Zin,c,modA,modF,Is,Isold,Isold2;
        rlgc,biasSI;
        width,height;
        v_TL,v_TLold;
        i_TL,i_TLold,i_TLold2;
        Bcoeff,Ecoeff,Fcoeff,coeff0fBC;
        Acoeff,Ccoeff,Dcoeff,Gcoeff;
        N_pts;
        IDX;
        
    end
    
    methods
        function obj = TL_solver3(params,rlgc)
            
        %Definition of constants
        obj.dx = params.dx;              %unit: mm
        obj.dt = params.dt;              %unit: picosecond
        obj.c = params.c;                %unit: mm/picosecond
        obj.IDX = params.IDX;
        obj.N_pts = params.N_pts;
               
        obj.rlgc = rlgc;
        obj.width = params.width;               %unit: mm
        obj.height = params.height;             %unit: mm
        
        obj.biasSI = params.bias/(obj.height*1e+3); %unit: kV/mm
        obj.modA = params.modA; %modulation amplitude unit: kV/mm
        obj.modF = params.modF;          %modulation frequency
        obj.Zin = params.Zin;            %unit: ohm
        obj.rlgc.R = 45*sqrt(obj.modF);  %unit: Ohm/mm --from paper W. Maineult: Microwave modulation of terahertz quantum cascade lasers: a transmission-line approach

        obj.Bcoeff = obj.dt/obj.dx/obj.rlgc.C*obj.width/(obj.height*1e+3);   % for calculation of v_TL
        obj.coeff0fBC = obj.Zin*obj.rlgc.C*obj.dx;
        
        %with consideration of resistance
        obj.Ecoeff = 2*obj.rlgc.L+obj.rlgc.R*obj.dt;
        obj.Fcoeff = 1/(obj.dx*(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R))*(obj.height*1e+3/obj.width);% *(V(2)-V(1))
        
        obj.Acoeff = -2/obj.dx^2+2*obj.rlgc.C*obj.rlgc.L/obj.dt^2;
        obj.Ccoeff = -obj.rlgc.C*obj.rlgc.L/obj.dt^2+obj.rlgc.C*obj.rlgc.R/2/obj.dt;
        obj.Dcoeff = obj.rlgc.C*obj.rlgc.L/obj.dt^2+obj.rlgc.C*obj.rlgc.R/2/obj.dt;
        
        %Initial conditions
        obj.v_TL = obj.biasSI*ones(obj.N_pts,1);    % V(0)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);            % i(1/2).... i(M+1/2)
        end
        

        
       function propagate3(obj,J_TL1,J_TL0,t)
  
        obj.v_TLold = obj.v_TL(1); % store the boundary voltage of previous step      
        obj.v_TL(1)= obj.v_TL(1)+obj.Bcoeff*(obj.Is-obj.i_TL(1)-(J_TL1(1)+J_TL0(1))*obj.dx/2); 
%         obj.v_TL(1) = obj.biasSI;
        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(1:end-2)-obj.i_TL(2:end-1)-(J_TL1(2:end-1)+J_TL0(2:end-1))*obj.dx/2);
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end-1)-obj.i_TL(end)-(J_TL1(end)+J_TL0(end))*obj.dx/2);

        
        
        obj.i_TLold = obj.i_TL; % store the boundary voltage of previous step 
        obj.i_TL(1) = obj.i_TL(1)-obj.dt/obj.dx/obj.rlgc.L*(obj.v_TL(2)-obj.v_TL(1))*(obj.height*1e+3)/obj.width;
     
        obj.i_TL(2:end-1) = (obj.Acoeff*obj.i_TL(2:end-1)...
                            +obj.Ccoeff*obj.i_TLold2(2:end-1)...
                            +1/obj.dx^2*(obj.i_TL(3:end)+obj.i_TL(1:end-2))...
                            +1/2/obj.dx*(J_TL1(3:end)-J_TL0(3:end)-J_TL1(2:end-1)+J_TL0(2:end-1)))/obj.Dcoeff;
        obj.i_TLold2 = obj.i_TLold;
        
        
        obj.Isold = obj.Is;
        obj.Is = obj.Isold2 + 1/obj.Zin/obj.width...
                 *(4*pi*obj.modF*obj.dt*obj.modA*cos(2*pi*obj.modF*t)...
                 -2*(obj.v_TL(1)-obj.v_TLold)*(obj.height*1e+3));
        obj.Isold2 = obj.Isold;                      
 
        end

        
          
    end
    
end

