classdef DM_MODEL_2_LVL_RWA_FP < handle
    % DM_MODEL_2_LVL_RWA_FP
    % Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        name;
        % upper laser level vars
        rho_u_DC; rho_u_DC_solver; rho_u_DC_t;
        rho_u_AC; rho_u_AC_solver; rho_u_AC_t;
        rho_u_DC_eq; % equilibrium values
        
        rho_l_DC; rho_l_DC_solver; rho_l_DC_t;
        rho_l_AC; rho_l_AC_solver; rho_l_AC_t;
        rho_l_DC_eq; % equilibrium values
        
        
        
        eta_ul_p; eta_ul_p_solver; eta_ul_p_t;
        eta_ul_m; eta_ul_m_solver; eta_ul_m_t;
        
        
        E_UL,E_LL, E0;
        
        T1,T2,DIFFUSION;
        
        IDX;
        zNORM;
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
        hbar;
        Gamma,zUL;
        Ncarriers;
        loss,nTHz;
        tch,lch; % characteristic time and length
    end
    
    properties(Constant)
        CAVITY = 'FP'
    end
    
    
    
    
    methods(Static)
        function params_gain = calc_threshold_inv(params_gain,params_abs,pump_strength)
            
            %The following section calculates the threshold inversion and sets the
            %correspoding steady state values for rho_u and rho_l
            
            %gain section carrier density (doping density) and trace normalization constant
            Ncarriers_g = params_gain.Ncarriers_cm*(100^3); % cm^-3 --> m^-3; carrier density
            %absorption section carrier density and trace normalization constant
            Ncarriers_a = params_abs.Ncarriers_cm*(100^3); % cm^-3 --> m^-3; carrier density
            % dipole moments
            mu_g = params_gain.zUL*1E-9*Constants('q0');
            mu_a = params_abs.zUL*1E-9*Constants('q0');
            % gain and absorber section cross sections
            
            % convert from eV to 2pi/ps (ang freq).
            hbar = Constants('hbar',{'time',params_gain.tch})/Constants('q0');
            E0 = params_gain.E0/hbar;
            
            sigma_g = params_gain.T_2*(1e-12)*params_gain.Gamma*E0*1E12*mu_g^2/ ...
                (Constants('eps0')*params_gain.nTHz*Constants('c')*Constants('hbar'));
            sigma_a = params_abs.T_2*(1e-12)*params_abs.Gamma*E0*1E12*mu_a^2/ ...
                (Constants('eps0')*params_abs.nTHz*Constants('c')*Constants('hbar'));
            %combined...
            Ga = sigma_a*Ncarriers_a*params_abs.L*1e-3;
            Gg = sigma_g*Ncarriers_g*params_gain.L*1e-3;
            powerloss = 2*params_gain.linear_loss*100;
            
            % threshold condition ! Note that here we assume tha the absorber is
            % a perfectly inverted two level system; i.e. rho_uu - rho_ll = -1 (at steady state)
            d_th =  powerloss/(sigma_g*Ncarriers_g)+Ga/Gg;
            params_gain.d_th = d_th;
            % either externally
            p = pump_strength;
            
            % how to set the inversion?!
            inv = p*d_th.* (p*d_th <=1 ) + 1*(p*d_th > 1);
            
            % steady state values for the upper level density and the lower level density;
            params_gain.rho_u_0 = (1+inv)/2;
            params_gain.rho_l_0 = (1-inv)/2;            
        end
    end
    
    methods
        function obj = DM_MODEL_2_LVL_RWA_FP(params)
            
            
            obj.name = params.name;
            %characteristic length and time
            obj.tch = params.tch;
            obj.lch = params.lch;
            
            %eigenenergies and dephasing times
            
            obj.hbar =  Constants('hbar',{'time',obj.tch})/Constants('q0');
            obj.E_UL = params.E_UL/obj.hbar; % in units of 2pi/time
            obj.E_LL = params.E_LL/obj.hbar; % in units of 2pi/time
            obj.E0 = params.E0/obj.hbar; % central freq in units 2pi/time;
            
            % upper level life time and dephasing time :)
            obj.T1 = params.T_1;
            obj.T2 = params.T_2;
            % diffusion time
%             obj.DIFFUSION =  (4*obj.E0/params.c^2*params.D*10^2/(1/params.tch));
            obj.DIFFUSION = 0;
            
            obj.rho_u_DC_eq = params.rho_u_0;
            obj.rho_l_DC_eq = params.rho_l_0;
            
            assert((obj.rho_u_DC_eq+obj.rho_l_DC_eq)==1,...
                'ASSERTION FAILURE! STEADY STATE VALUE IS WRONG')
            
            % nomalization constant params
            obj.Gamma = params.Gamma;  % overlap factor
            obj.zUL = params.zUL ; % dipole element in units of nm
            obj.nTHz = params.nTHz;
            
            obj.N_pts = params.N_pts;
            % initialize the variables
            obj.rho_u_DC = obj.rho_u_DC_eq*ones(obj.N_pts,1);
            obj.rho_l_DC = obj.rho_l_DC_eq*ones(obj.N_pts,1);
            obj.rho_u_AC = zeros(obj.N_pts,1);
            obj.rho_l_AC = zeros(obj.N_pts,1);
            
            
            %some very small correlation
            obj.eta_ul_p = 1e-10*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.eta_ul_m = 1e-10*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            
            obj.eta_ul_p_t = 0*obj.eta_ul_p;
            obj.eta_ul_m_t = 0*obj.eta_ul_m;
            
            
            
            obj.rho_u_DC_solver = MS(5,obj.N_pts,[],obj.rho_u_DC);
            obj.rho_u_AC_solver = MS(5,obj.N_pts,[],obj.rho_u_AC);
            
            obj.rho_l_DC_solver = MS(5,obj.N_pts,[],obj.rho_l_DC);
            obj.rho_l_AC_solver = MS(5,obj.N_pts,[],obj.rho_l_AC);
            
            obj.eta_ul_p_solver = MS(5,obj.N_pts,[],obj.eta_ul_p);
            obj.eta_ul_m_solver = MS(5,obj.N_pts,[],obj.eta_ul_m);
            
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
            
        end
        function [] = propagate(obj,F,B,dt) % timestep the density matrix
            %%%
            % F: forward propagating component of the field
            % B: backward propagating component of the field
            % dt: time-step
            %%%
            
            norm_ftr = obj.NORM_FACTOR_DM;
            
            % light matter interaction ->
            lmInteraction = norm_ftr*(conj(F).*obj.eta_ul_p+conj(B).*obj.eta_ul_m);
            
            %populations
            %%% upper laser level
            obj.rho_u_DC_t = 1i/2*(lmInteraction-conj(lmInteraction)) - ...
                (obj.rho_u_DC-obj.rho_u_DC_eq)/obj.T1;
            obj.rho_u_DC_solver.make_step(obj.rho_u_DC_t,dt);
            
            obj.rho_u_AC_t = 1i/2*norm_ftr*(conj(B).*obj.eta_ul_p-F.*conj(obj.eta_ul_m)) -...
                (1/obj.T1+obj.DIFFUSION)*obj.rho_u_AC;
            obj.rho_u_AC_solver.make_step(obj.rho_u_AC_t,dt);
            
            %%% lower laser level
            obj.rho_l_DC_t = -1i/2*(lmInteraction-conj(lmInteraction)) - ...
                (obj.rho_l_DC-obj.rho_l_DC_eq)/obj.T1;
            obj.rho_l_DC_solver.make_step(obj.rho_l_DC_t,dt);
            
            obj.rho_l_AC_t = -1i/2*norm_ftr*(conj(B).*obj.eta_ul_p-F.*conj(obj.eta_ul_m)) - ...
                (1/obj.T1+obj.DIFFUSION)*obj.rho_l_AC;
            obj.rho_l_AC_solver.make_step(obj.rho_l_AC_t,dt);
            
            % inversion
            D0 = obj.rho_u_DC-obj.rho_l_DC;
            % population grating
            D2 = obj.rho_u_AC-obj.rho_l_AC;
            
            %coherences
            obj.eta_ul_p_t = (1i*(obj.E0-(obj.E_UL-obj.E_LL))-1/obj.T2).*obj.eta_ul_p + ...
                1i/2*norm_ftr.*(F.*D0+B.*D2);
            obj.eta_ul_p_solver.make_step(obj.eta_ul_p_t,dt);
            
            obj.eta_ul_m_t = (1i*(obj.E0-(obj.E_UL-obj.E_LL))-1/obj.T2).*obj.eta_ul_m + ...
                1i/2*norm_ftr.*(B.*D0+F.*conj(D2));
            obj.eta_ul_m_solver.make_step(obj.eta_ul_m_t,dt);
            
        end
        
        function [] = update_state(obj)
            obj.rho_u_DC = obj.rho_u_DC_solver.get_latest_solution();
            obj.rho_u_AC = obj.rho_u_AC_solver.get_latest_solution();
            
            obj.rho_l_DC = obj.rho_l_DC_solver.get_latest_solution();
            obj.rho_l_AC = obj.rho_l_AC_solver.get_latest_solution();
            
            obj.eta_ul_p = obj.eta_ul_p_solver.get_latest_solution();
            obj.eta_ul_m = obj.eta_ul_m_solver.get_latest_solution();
        end
        
        function [P,P_t,M,M_t,LOSSES] = get_polarization_and_losses(obj)
            % additionally the polarizaiton should be multiplied by -1i*c/n
            P = obj.NORM_FACTOR_FIELD*obj.eta_ul_p;
            P_t = obj.NORM_FACTOR_FIELD*obj.eta_ul_p_t;
            
            M = obj.NORM_FACTOR_FIELD*obj.eta_ul_m;
            M_t = obj.NORM_FACTOR_FIELD*obj.eta_ul_m_t;
            
            LOSSES = obj.loss;
        end
        
        function [gn] = calculate_optical_gain(obj,nthz)
            
            avg_inv = mean(obj.rho_u_DC-obj.rho_l_DC);
            
            gn = (obj.E0*obj.T2/(Constants('hbar')*Constants('eps0')*...
                Constants('c')*nthz))*obj.Ncarriers*...
                (obj.zUL*1e-9*Constants('q0'))^2*avg_inv;
            
        end
        
         function [gn] = calculate_positional_gain(obj,nthz)
            gn = (obj.E0*obj.T2/(Constants('hbar')*Constants('eps0')*...
                Constants('c')*nthz))*obj.Ncarriers*...
                (obj.zUL*1e-9*Constants('q0'))^2*(obj.rho_u_DC-obj.rho_l_DC);
            
        end
        
        
    end
    
end

