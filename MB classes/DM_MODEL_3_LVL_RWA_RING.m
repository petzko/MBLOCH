classdef DM_MODEL_3_LVL_RWA_RING < handle
    %DM_MODEL_AC_M_LVL_RWA_RING Summary of this class goes here
    %   Detailed explanation goes here
    properties
        name;
        
        % upper laser level state
        rho_u; rho_u_solver;
        rho_u_t;
        
        % injector state (RT-coupled  to the excited state)
        rho_i; rho_i_solver;
        rho_i_t;
        
        % lower laser level state
        rho_l; rho_l_solver;
        rho_l_t;
        
        
        % coherences
        eta_ul; eta_il; rho_iu;
        eta_ul_t; eta_il_t; rho_iu_t;
        eta_ul_solver; eta_il_solver; rho_iu_solver;
        
        % Hamiltonian in units (1/time)
       % Energies in units (1/time)
        E_UL, E_LL, E_IL, 
        % anticrossing frequency Omega and field central frequency E0
        Omega; E0;hbar;
        
        % scattering rates smatrix
        W; 
        % inversion recovery dephasing times; 
        T1;T2; 
        % 1/lifetime array!
        TOT_RTS; DEPH; IDX;
        
        % Indices of the ULL, INJ and the LLL  state
        ULL, INJ, LLL, DEPOP;
        
        % indices of the different polarization terms int he deph array
        idx_iu, idx_ul,idx_il;                
        NORM_FACTOR_FIELD; % factor that I have to multiply the eta_ul with
        % in the field propagation equations;
        NORM_FACTOR_DM; % factor that I have to multiply the electric field in
        
        % the density matrix equations.
        % If the global electric field is
        % normalized w.r.t. some dipole element zNORM, i.e. E_norm =
        % zNORM*E/hbar ->
        % NORM_FACTOR_FIELD = TRACE_RHO* zNORM/zUL
        % and NORM_FACTOR_DM =  zUL/zNORM;
        
        TRACE_RHO; % =  Ncarriers*Gamma*E0*(qzUL^2)/eps0*nthz*c*hbar
        
        N_pts;NLVLS;
        Gamma,zUL,zNORM;
        Ncarriers;
        loss,nTHz;
        tch,lch; % characteristic time and length
    end
    
    methods

        function obj  = DM_MODEL_3_LVL_RWA_RING(params)
            
            obj.name = params.name;
            %characteristic length and time
            obj.tch = params.tch;
            obj.lch = params.lch;
            
            obj.hbar =  Constants('hbar',{'time',obj.tch})/Constants('q0');
            obj.E_UL = params.E_UL/obj.hbar; % in units of 2pi/time
            obj.E_LL = params.E_LL/obj.hbar; % in units of 2pi/time
            obj.E_IL = params.E_IL/obj.hbar; % in units of 2pi/time
            obj.Omega = params.AC_energy/obj.hbar;  % in units of 2pi/time
            obj.E0 = params.E0/obj.hbar; % central freq in units 2pi/time;
            
            
            % nomalization constant params
            obj.Gamma = params.Gamma;  % overlap factor
            obj.zUL = params.zUL ; % dipole element in units of nm
            obj.nTHz = params.nTHz;
            
            obj.N_pts = params.N_pts;
            obj.NLVLS = 4; 
            
            obj.ULL = params.ULL;
            obj.INJ = params.INJ;
            obj.LLL = params.LLL;
            obj.DEPOP = params.DEPOP;
     
     
            obj.loss = params.linear_loss*100/(1/obj.lch)*ones(obj.N_pts,1);
            obj.IDX = params.IDX;
            
            % Carrier density in 1/m^3!
            obj.Ncarriers = params.Ncarriers_cm*(100^3);
            obj.TRACE_RHO = ((obj.E0*1e12*obj.Ncarriers*obj.Gamma*...
                (obj.zUL*1E-9*Constants('q0'))^2)...
                /(Constants('eps0')*obj.nTHz*Constants('c')*Constants('hbar')))...
                /(1/(obj.lch*obj.tch));

            obj.zNORM = params.zNORM;
            obj.NORM_FACTOR_FIELD = obj.TRACE_RHO*obj.zNORM/obj.zUL;
            obj.NORM_FACTOR_DM = obj.zUL/obj.zNORM;
              
            % initialize the variables
            obj.rho_u = 1/3*ones(obj.N_pts,1);
            obj.rho_i = 1/3*ones(obj.N_pts,1);
            obj.rho_l =  1-obj.rho_u-obj.rho_i;
            
            % some very small correlation
            obj.eta_ul = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.rho_iu = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.eta_il = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            
            obj.eta_ul_t = 0*obj.eta_ul;
            obj.eta_il_t = 0*obj.eta_il;
            obj.rho_iu_t = 0*obj.rho_iu;
            
            obj.rho_i_solver = MS(5,obj.N_pts,[],obj.rho_i);
            obj.rho_u_solver = MS(5,obj.N_pts,[],obj.rho_u);
            obj.rho_l_solver = MS(5,obj.N_pts,[],obj.rho_l);
            
            obj.rho_iu_solver = MS(5,obj.N_pts,[],obj.rho_iu);
            obj.eta_ul_solver = MS(5,obj.N_pts,[],obj.eta_ul);
            obj.eta_il_solver = MS(5,obj.N_pts,[],obj.eta_il);
             
                            
            obj.DEPH = zeros(3,1);
            obj.W = zeros(4,4);
            obj.TOT_RTS = zeros(3,1);
            
            obj.idx_iu =1;
            obj.idx_ul =2; 
            obj.idx_il = 3;
            
          
            for i_ = 1:obj.NLVLS
                for j_ = 1:obj.NLVLS
                    if i_ == j_
                        continue
                    end
                    obj.W(i_,j_) = params.W(j_+(i_-1)*obj.NLVLS);
                end
            end
            
            obj.TOT_RTS(obj.INJ) = 0*obj.TOT_RTS(obj.INJ);
            for el = setdiff(1:obj.NLVLS,[obj.INJ,obj.DEPOP])
                obj.TOT_RTS(obj.INJ) = obj.TOT_RTS(obj.INJ)+obj.W(obj.INJ,el);
            end
            for lvl = setdiff(1:obj.NLVLS,[obj.INJ,obj.DEPOP])
                obj.TOT_RTS(lvl) = 0*obj.TOT_RTS(lvl);
                for el = 1:obj.NLVLS
                    if el ~=lvl
                        obj.TOT_RTS(lvl) = obj.TOT_RTS(lvl) +  obj.W(lvl,el);
                    end
                end
            end
            
            % Added Pure Dephasing
            obj.DEPH(obj.idx_iu) = 1/2*(obj.TOT_RTS(obj.INJ)+ ...
                obj.TOT_RTS(obj.ULL)) +1/params.Tdeph_iu;
            obj.DEPH(obj.idx_ul) = 1/2*(obj.TOT_RTS(obj.ULL)+...
                obj.TOT_RTS(obj.LLL))+1/params.Tdeph_ul;
            obj.DEPH(obj.idx_il) = 1/2*(obj.TOT_RTS(obj.INJ)+...
                obj.TOT_RTS(obj.LLL)) +1/params.Tdeph_il;
            obj.T2 = obj.DEPH(obj.idx_ul); 
            obj.T1 = 1./obj.W(obj.ULL,obj.LLL); % approximately ! 
        end
        
        function [] = propagate(obj,FIELD,dt)
            % POPULATIONS
            ULL = obj.ULL; INJ = obj.INJ;  LLL = obj.LLL; DEPOP = obj.DEPOP;
            norm_ftr = obj.NORM_FACTOR_DM;
          
            % rho_i
            obj.rho_i_t = 1i*obj.Omega.*(obj.rho_iu- ...
                conj(obj.rho_iu))+obj.W(ULL,INJ).*obj.rho_u +...
                obj.W(LLL,INJ)*obj.rho_l - obj.TOT_RTS(INJ)*obj.rho_i;
            obj.rho_i_solver.make_step(obj.rho_i_t,dt);
            
            % rho_u
            lmInteraction = norm_ftr.*conj(FIELD).*obj.eta_ul;
            obj.rho_u_t = -1i*obj.Omega.*(obj.rho_iu- ...
                conj(obj.rho_iu)) + 1i/2.*(lmInteraction-conj(lmInteraction)) +...
                obj.rho_i*(obj.W(INJ,ULL)) +...
                obj.rho_l.*obj.W(LLL,ULL) - obj.TOT_RTS(ULL)*obj.rho_u;
            obj.rho_u_solver.make_step(obj.rho_u_t,dt);
            
            % rho_l
            obj.rho_l_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + ...
                obj.rho_i.*(obj.W(INJ,LLL)) +...
                obj.rho_u.*obj.W(ULL,LLL) - obj.TOT_RTS(LLL)*obj.rho_l;
            obj.rho_l_solver.make_step(obj.rho_l_t,dt);
                                 
            %%% coherences! rho_iu, eta_ul and eta_il
            % rho_iu
            obj.rho_iu_t = (-1i*(obj.E_IL-obj.E_UL) - obj.DEPH(obj.idx_iu)).*obj.rho_iu + ...
                1i*obj.Omega.*(obj.rho_i - obj.rho_u) + ...
                1i/2*norm_ftr.*conj(FIELD).*obj.eta_il;
            obj.rho_iu_solver.make_step(obj.rho_iu_t,dt);
            
            % eta_ul
            obj.eta_ul_t = (1i*(obj.E0-(obj.E_UL-obj.E_LL)) - obj.DEPH(obj.idx_ul)).*obj.eta_ul + ...
                1i/2*norm_ftr.*FIELD.*(obj.rho_u-obj.rho_l)-1i*obj.Omega.*obj.eta_il;
            obj.eta_ul_solver.make_step(obj.eta_ul_t,dt);
            
            % eta_il
            obj.eta_il_t = (1i*(obj.E0-(obj.E_IL-obj.E_LL)) - obj.DEPH(obj.idx_il))./obj.eta_il + ...
                1i/2*norm_ftr.*FIELD.*obj.rho_iu -1i*obj.Omega.*obj.eta_ul;
            obj.eta_il_solver.make_step(obj.eta_il_t,dt);
        end
        
        function tr = get_avg_trace(obj)
            tr = obj.rho_i+obj.rho_u+obj.rho_l; 
            tr = mean(tr);        
        end
        
        function [] = update_state(obj)
            obj.rho_l = obj.rho_l_solver.get_latest_solution();
            obj.rho_u = obj.rho_u_solver.get_latest_solution();
            obj.rho_i = obj.rho_i_solver.get_latest_solution();
 
            obj.rho_iu = obj.rho_iu_solver.get_latest_solution();
            obj.eta_ul = obj.eta_ul_solver.get_latest_solution();
            obj.eta_il = obj.eta_il_solver.get_latest_solution();
        end
        
        function [P,P_t,LOSSES] = get_polarization_and_losses(obj)
            % additionally the polarizaiton should be multiplied by -1i*c/n
            P = obj.NORM_FACTOR_FIELD*obj.eta_ul;
            P_t = obj.NORM_FACTOR_FIELD*obj.eta_ul_t;
            LOSSES = obj.loss;
        end
        
    end
    
end