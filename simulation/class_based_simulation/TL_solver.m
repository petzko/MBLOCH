classdef TL_solver < handle
    %TL_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dx,dt,Zin,Vs;
        rlgc;
        width,height;
        v_TL,i_TL,J_TL;
        N_pts;
        IDX;
        
    end
    
    methods
        function obj = TL_solver(params,rlgc)
            
        %Definition of constants
        obj.dx = params.dx;        %unit: mm
        obj.dt = params.dt*1e-12;        %unit: ps
        obj.Zin = params.Zin;      %unit: ohm
        obj.Vs = params.Vs;        %unit: V
        obj.IDX = params.IDX;
        obj.N_pts = params.N_pts;
        
        obj.rlgc = rlgc;

        obj.width = params.width;     %unit: mm
        obj.height = params.height;     %unit: mm
     
        %Initial conditions
        %At t=0, the voltage on lines were set to be equal with half of the source,
        %and no current flows to lines. Due to tunneling of electrons (J) which is 
        %relevant to position, the voltage as well as current along the lines 
        %will be changed until equilibrium.
        obj.v_TL = obj.Vs/2*ones(obj.N_pts,1);    % V(0)    V(1)....   V(M)
        obj.i_TL = 0*ones(obj.N_pts,1);    % I(1/2)  I(3/2).... I(M+1/2)
        end
        
        function propagate(obj,J_TL)
            
        obj.J_TL = J_TL;

        obj.v_TL(2:end-1) = obj.v_TL(2:end-1)-obj.dt/obj.rlgc.C/obj.dx*(obj.i_TL(2:end-1)-...
                            obj.i_TL(1:end-2)+obj.width*obj.dx*J_TL(2:end-1));
        obj.i_TL(1:end-1) = obj.i_TL(1:end-1)-obj.dt/obj.rlgc.L/obj.dx*(obj.v_TL(2:end)...
                            -obj.v_TL(1:end-1));

        end
       
        function set_boundary(obj,v2old,v2new,i2new,width2,J_2TL,rlgc2)
            
%             obj.v_TL(1) = 2*(obj.Vs-(v2old+v2new)/2-obj.Zin*(i2new+width2*obj.dx*J_2TL(1)+...
%                             rlgc2.C*obj.dx/obj.dt/2*(v2new-v2old)))-obj.v_TL(1);
            obj.v_TL(1) = 1/(rlgc2.C*obj.dx/obj.dt/2+1/2/obj.Zin)*...
                    ((rlgc2.C*obj.dx/obj.dt/2-1/2/obj.Zin)*v2new-i2new-width2*obj.dx/2*...
                    J_2TL(1)+(obj.Vs-(v2new+v2old)/2)/obj.Zin);
            obj.i_TL(end) = -obj.i_TL(end-1);
            obj.v_TL(end) = obj.v_TL(end)-obj.dt/obj.rlgc.C/obj.dx*(-2*obj.i_TL(end-1)...
                        +obj.width*obj.dx*J_2TL(end));
        end
        
        
        function v = get_bias(obj)
            v = obj.v_TL*1e-3/obj.height;
        end
          
    end
    
end

