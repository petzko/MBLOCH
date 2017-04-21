classdef TL_solver2 < handle
    %TL_SOLVER Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        dx,dt,Zin,c,modA,modF,V_RFnew,V_RFold;
        rlgc,bias;
        width,height;
        v_TL,v_TLold,v_TLold2,Mirror_v;
        i_TL,i_TLold,i_TLold2,Mirror_i;
        Bcoeff,,Ecoeff,Fcoeff,coeff0fBC;
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
        obj.coeff0fBC = obj.Zin*obj.rlgc.C*obj.dx/obj.dt;
        
        %with consideration of resistance
        obj.Ecoeff = obj.rlgc.L*obj.dx/obj.dt;
        obj.Fcoeff = -1/(obj.dx*(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R))*(obj.height*1e+3/obj.width);% *(V(2)-V(1))
        
        %Initial conditions
        obj.v_TL = obj.bias*ones(obj.N_pts,1);    % V(0)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);            % i(1/2).... i(M+1/2)
        obj.v_TLold2 = obj.v_TL; % previous voltage at node 0 before initial (t<0)
        obj.V_RFnew = 0;
        obj.V_RFold = 0;
        end
        

        
       function propagate2(obj,J_TL1,J_TL0,t)
       
        obj.v_TLold = obj.v_TL; % store the boundary voltage of previous step      
        obj.v_TL(1)=(-obj.coeff0fBC*obj.v_TL(1)+(1+obj.coeff0fBC)/2*obj.v_TLold2(1)+obj.V_RFnew-obj.V_RFold...
                      + (obj.i_TL(1)-obj.i_TLold(1))*obj.width*obj.Zin/(obj.height*1e3))/(1/2-obj.coeff0fBC/2);          
        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(2:end-1)-obj.i_TL(1:end-2)+obj.dx*((J_TL1(2:end-1)+J_TL0(2:end-1))/2));
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end)-obj.i_TL(end-1)+(J_TL1(end)+J_TL0(end))*obj.dx/4)*2;
        obj.v_TLold2 = obj.v_TLold;
        
        
        obj.i_TLold = obj.i_TL; % store the boundary voltage of previous step      
        obj.i_TL(1:end-1) = (2*obj.Ecoeff*obj.i_TL(1:end-1)+(-obj.Ecoeff+obj.rlgc.R*obj.dx/2)*obj.i_TLold2(1:end-1)-...
                            obj.height*1e+3/obj.width*(obj.v_TL(2:end)-obj.v_TL(1:end-1)-obj.v_TLold(2:end)+obj.v_TLold(1:end-1)))/...
                            (obj.Ecoeff+obj.rlgc.R*obj.dx/2);
        obj.i_TLold2 = obj.i_TLold;
 
        
        obj.V_RFnew = obj.modA*sin(2*pi*obj.modF*(t-obj.dt/2))/(obj.height*1e+3);    
        obj.V_RFold = obj.modA*sin(2*pi*obj.modF*(t-obj.dt*1.5))/(obj.height*1e+3);
        end

        
          
    end
    
end

