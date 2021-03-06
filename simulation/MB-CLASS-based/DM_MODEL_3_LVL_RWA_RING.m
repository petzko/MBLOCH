classdef DM_MODEL_3_LVL_RWA_RING < handle
    %DM_MODEL_3_LVL_RWA_RING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        
        % excited state
        rho_ex; rho_ex_solver;
        rho_ex_t;
        
        % spin state (RT-coupled  to the excited state)
        rho_spin; rho_spin_solver;
        rho_spin_t;
        
        % ground state
        rho_grnd; rho_grnd_solver;
        rho_grnd_t;
        
        % coherences
        eta_eg; eta_sg; rho_se;
        eta_eg_t; eta_sg_t; rho_se_t;
        eta_eg_solver; eta_sg_solver;rho_se_solver;
        
        % Hamiltonian in units (1/time)
        H;
        E0;
        
        
        % Indices of the Excited, Spin and the Ground state
        Ex_, Spin_, Grnd_;
        % indices of the different polarization terms int he deph array
        SE, EG,SG;
        
        % scattering rate smatrix
        W;
        % 1/lifetime array!
        TOT_RTS;
        DEPH;
        
        IDX;
        
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
        
        N_pts;
        Gamma,zUL;
        Ncarriers;
        loss,nTHz;
        tch,lch; % characteristic time and length
    end
    
    properties (Constant)
        CAVITY = 'RING'
    end
    
    methods
        
        function obj  = DM_MODEL_3_LVL_RWA_RING(params)
            
            
            assert(strcmp(params.cavity,obj.CAVITY),'DM MODEL ERROR! CAVITY TYPE MISMATCH.');
            
            obj.name = params.name;
            %characteristic length and time
            obj.tch = params.tch;
            obj.lch = params.lch;
            
            %eigenenergies and dephasing times
            
            hbar =  Constants('hbar',{'time',obj.tch})/Constants('q0');
            obj.H = params.H/hbar; % in units of 2pi/time
            obj.E0 = params.E0/hbar;
            
            obj.W = params.W;
            
            % nomalization constant params
            obj.Gamma = params.Gamma;  % overlap factor
            obj.zUL = params.zUL ; % dipole element in units of nm
            obj.nTHz = params.nTHz;
            
            obj.N_pts = params.N_pts;
            % initialize the variables
            obj.rho_ex = 1/3*ones(obj.N_pts,1);
            obj.rho_spin = 1/3*ones(obj.N_pts,1);
            obj.rho_grnd =  1 - obj.rho_ex-obj.rho_spin;
            
            % some very small correlation
            obj.eta_eg = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.rho_se = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.eta_sg = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            
            obj.eta_eg_t = 0*obj.eta_eg;
            obj.eta_sg_t = 0*obj.eta_sg;
            obj.rho_se_t = 0*obj.rho_se;
            
            obj.rho_spin_solver = MS(5,obj.N_pts,[],obj.rho_spin);
            obj.rho_ex_solver = MS(5,obj.N_pts,[],obj.rho_ex);
            obj.rho_grnd_solver = MS(5,obj.N_pts,[],obj.rho_grnd);
            
            obj.rho_se_solver = MS(5,obj.N_pts,[],obj.rho_se);
            obj.eta_eg_solver = MS(5,obj.N_pts,[],obj.eta_eg);
            obj.eta_sg_solver = MS(5,obj.N_pts,[],obj.eta_sg);
            
            obj.Ex_ = params.Ex_;
            obj.Spin_ = params.Spin_;
            obj.Grnd_ = params.Grnd_;
            
            
            obj.SE = 1; obj.EG = 2; obj.SG = 3;
            
            obj.TOT_RTS = zeros(3,1);
            obj.TOT_RTS(obj.Spin_) = obj.W(obj.Spin_,obj.Grnd_) + ...
                obj.W(obj.Spin_,obj.Ex_);
            obj.TOT_RTS(obj.Ex_) = obj.W(obj.Ex_,obj.Grnd_) + ...
                obj.W(obj.Ex_,obj.Spin_);
            obj.TOT_RTS(obj.Grnd_) = obj.W(obj.Grnd_,obj.Spin_) + ...
                obj.W(obj.Grnd_,obj.Ex_);
            
            obj.DEPH = zeros(3,1);
            obj.DEPH(obj.SE) = 1/2*(obj.TOT_RTS(obj.Spin_)+obj.TOT_RTS(obj.Ex_));
            obj.DEPH(obj.EG) = 1/2*(obj.TOT_RTS(obj.Ex_)+obj.TOT_RTS(obj.Grnd_));
            obj.DEPH(obj.SG) = 1/2*(obj.TOT_RTS(obj.Spin_)+obj.TOT_RTS(obj.Grnd_));
            
            
            obj.loss = params.linear_loss*100/(1/obj.lch)*ones(obj.N_pts,1);
            obj.IDX = params.IDX;
            
            % Carrier density in 1/m^3!
            obj.Ncarriers = params.Ncarriers_cm*(100^3);
            obj.TRACE_RHO = ((obj.E0*1e12*obj.Ncarriers*obj.Gamma*...
                (obj.zUL*1E-9*Constants('q0'))^2)...
                /(Constants('eps0')*obj.nTHz*Constants('c')*Constants('hbar')))...
                /(1/(obj.lch*obj.tch));
            obj.NORM_FACTOR_FIELD = obj.TRACE_RHO*params.zNORM/obj.zUL;
            obj.NORM_FACTOR_DM = obj.zUL/params.zNORM;
            
        end
        
        function [] = propagate(obj,FIELD,dt)
            % POPULATIONS
            E_idx = obj.Ex_; S_idx = obj.Spin_; G_idx = obj.Grnd_;
            norm_ftr = obj.NORM_FACTOR_DM;
            
            obj.rho_spin_t = 1i*obj.H(S_idx,E_idx).*(obj.rho_se- ...
                conj(obj.rho_se))+obj.W(E_idx,S_idx).*obj.rho_ex +...
                obj.W(G_idx,S_idx)*obj.rho_grnd - obj.TOT_RTS(S_idx)*obj.rho_spin;
            obj.rho_spin_solver.make_step(obj.rho_spin_t,dt);
            
            lmInteraction = norm_ftr.*conj(FIELD).*obj.eta_eg;
            obj.rho_ex_t = -1i*obj.H(S_idx,E_idx).*(obj.rho_se- ...
                conj(obj.rho_se)) + 1i/2.*(lmInteraction-conj(lmInteraction)) +...
                obj.rho_spin*(obj.W(S_idx,E_idx)) +...
                obj.rho_grnd.*obj.W(G_idx,E_idx) - obj.TOT_RTS(E_idx)*obj.rho_ex;
            obj.rho_ex_solver.make_step(obj.rho_ex_t,dt);
            
            obj.rho_grnd_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + ...
                obj.rho_spin.*(obj.W(S_idx,G_idx)) +...
                obj.rho_ex.*obj.W(E_idx,G_idx) - obj.TOT_RTS(G_idx)*obj.rho_grnd;
            obj.rho_grnd_solver.make_step(obj.rho_grnd_t,dt);
            
            %%% coherences! rho_se, eta_eg and eta_sg
            % rho_se
            obj.rho_se_t = -(1i*(obj.H(S_idx,S_idx)-obj.H(E_idx,E_idx)) + ...
                obj.DEPH(obj.SE)).*obj.rho_se + ...
                1i*obj.H(S_idx,E_idx).*(obj.rho_spin - obj.rho_ex) + ...
                1i/2*norm_ftr.*conj(FIELD).*obj.eta_sg;
            obj.rho_se_solver.make_step(obj.rho_se_t,dt);
            
            % eta_eg
            obj.eta_eg_t = -(1i*(obj.H(E_idx,E_idx)-obj.H(G_idx,G_idx)-obj.E0) + ...
                obj.DEPH(obj.EG)).*obj.eta_eg + ...
                1i/2*norm_ftr.*FIELD.*(obj.rho_ex-obj.rho_grnd) - ...
                1i*obj.H(S_idx,E_idx).*obj.eta_sg;
            obj.eta_eg_solver.make_step(obj.eta_eg_t,dt);
            
            % eta_sg
            obj.eta_sg_t = -(1i*(obj.H(S_idx,S_idx)-obj.H(G_idx,G_idx)-obj.E0)+ ...
                obj.DEPH(obj.SG)).*obj.eta_sg + 1i/2*norm_ftr.*FIELD.*obj.rho_se -...
                1i*obj.H(S_idx,E_idx).*obj.eta_eg;
            obj.eta_sg_solver.make_step(obj.eta_sg_t,dt);
        end
        
        function [] = update_state(obj)
            obj.rho_grnd = obj.rho_grnd_solver.get_latest_solution();
            obj.rho_ex = obj.rho_ex_solver.get_latest_solution();
            obj.rho_spin = obj.rho_spin_solver.get_latest_solution();
            
            obj.rho_se = obj.rho_se_solver.get_latest_solution();
            obj.eta_eg = obj.eta_eg_solver.get_latest_solution();
            obj.eta_sg = obj.eta_sg_solver.get_latest_solution();
        end
        
        
        function [P,P_t,LOSSES] = get_polarization_and_losses(obj)
            % additionally the polarizaiton should be multiplied by -1i*c/n
            P = obj.NORM_FACTOR_FIELD*obj.eta_eg;
            P_t = obj.NORM_FACTOR_FIELD*obj.eta_eg_t;
            LOSSES = obj.loss;
        end
        
    end
    
end