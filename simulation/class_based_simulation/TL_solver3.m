classdef TL_solver3 < handle
    %TL_SOLVER Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        dx,dt,c,modA,modF,Vs;
        rlgc,bias,wire;
        width,height;
        v_TL,v_TLold,v_TLold2;
        i_TL,i_TLold;
        Bcoeff,,Ecoeff,Fcoeff,coeff0fBC;
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
        obj.wire.C = 27/0.3048*1e-3;                        %unit: pF    pF/ft /0.3048 m/ft
        obj.wire.L = 6.7498e4/0.3048*1e-3;                  %unit: pH
        
        obj.bias = params.bias/(obj.height*1e+3); %unit: kV/mm
        obj.modA = params.modA/(obj.height*1e+3); %modulation amplitude unit: kV/mm
        obj.modF = params.modF;          %modulation frequency

        obj.rlgc.R = 0;%45*sqrt(obj.modF);  %unit: Ohm/mm --from paper W. Maineult: Microwave modulation of terahertz quantum cascade lasers: a transmission-line approach

        obj.Bcoeff = -obj.width/(obj.height*1e+3)/obj.c/obj.rlgc.C;   % for calculation of v_TL
        obj.coeff0fBC = obj.wire.L*obj.wire.C/(obj.dt)^2;
        
        %with consideration of resistance
        obj.Ecoeff = (obj.rlgc.L/obj.dt-1/2*obj.rlgc.R)/(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R);
        obj.Fcoeff = -1/(obj.dx*(obj.rlgc.L/obj.dt+1/2*obj.rlgc.R))*(obj.height*1e+3/obj.width);% *(V(2)-V(1))
        
        %Initial conditions
        obj.v_TL = obj.bias*ones(obj.N_pts,1);    % V(0)....   V(M)       
        obj.i_TL = zeros(obj.N_pts,1);            % i(1/2).... i(M+1/2)
        obj.v_TLold2 = obj.v_TL(1); % previous voltage at node 0 before initial (t<0)
        end
        

        
       function propagate3(obj,J_TL1,J_TL0,t)
           
        obj.Vs = obj.bias+obj.modA*sin(2*pi*obj.modF*(t-obj.dt))/(obj.height*1e+3);       
        obj.v_TLold = obj.v_TL(1); % store the boundary voltage of previous step
        
        obj.v_TL(1)= (2-1/obj.coeff0fBC)*obj.v_TL(1)-obj.v_TLold2+1/obj.coeff0fBC*...
                    (obj.Vs-obj.wire.L/obj.dt*(obj.i_TL(1)-obj.i_TLold)*obj.width/(obj.height*1e+3));
        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)+obj.Bcoeff*(obj.i_TL(2:end-1)-obj.i_TL(1:end-2)+(J_TL1(2:end-1)+J_TL0(2:end-1))*obj.dx/2);
        obj.v_TL(end) = obj.v_TL(end)+obj.Bcoeff*(obj.i_TL(end)-obj.i_TL(end-1)+(J_TL1(end)+J_TL0(end))*obj.dx/2)*2;
        obj.v_TLold2 = obj.v_TLold;
        
        obj.i_TLold = obj.i_TL(1); % store the boundary voltage of previous step      
        obj.i_TL(1:end-1) = obj.Ecoeff*obj.i_TL(1:end-1)+obj.Fcoeff*(obj.v_TL(2:end)-obj.v_TL(1:end-1));

        end

        
          
    end
    
end

