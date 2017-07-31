classdef DM_MODEL_MLVL_RWA_RING < handle
    properties
        name;
        
        % injector state (RT-coupled  to the excited state)
        rho_i;  rho_i_solver;  rho_i_t;
        
        % excited state
        rho_u;  rho_u_solver;   rho_u_t;
        
        % ground state
        rho_l; rho_l_solver;  rho_l_t;
        
        % rest of the levels~
        rRES; rRES_t; rRES_solver;
        % coherences backward and forward components.
        eta_ul;
        
        eta_ul_t; 
        
        
        eta_ul_solver; 
        
        
        % Energies in units (1/time)
        E_UL;
        dE_UL;
        % central frequency
        E0;
        % period length
        Lp;
        % total section length;
        L;
        
        % Indices of the Excited, Spin and the Ground state as well as the
        % depopulation level
        ULL, INJ, LLL, DEPOP,IGNORELEVEL;
        
        % total number of levels
        NLVLS; global_idx_rest;        N_rest;
        steady_state;
        
        % indices of the different polarization terms int he deph array
        i_ul
        
        % scattering rates matrix/reduced S matrix/ pump vector;
        W; S, b;
        
        d_th; % threshold inversion!
        
        % 1/lifetime array!
        TOT_RTS;
        DEPH;
        Tdeph_ul;
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
        Gamma,zUL,zNORM;
        Ncarriers;
        loss,nTHz;
        tch,lch; % characteristic time and length
        
        
        
    end
    
    
    methods
        
        function obj  = DM_MODEL_MLVL_RWA_RING(params)
            
            obj.name = params.name;
            %characteristic length and time
            obj.tch = params.tch;
            obj.lch = params.lch;
            
            % ALLlevels
            obj.NLVLS = params.NLVLS;
            obj.Lp = params.Lp;
            obj.L = params.L;
            % nomalization constant params
            obj.Gamma = params.Gamma;  % overlap factor
            
            
            obj.nTHz = params.nTHz;
            obj.N_pts = params.N_pts;
            obj.ULL = params.ULL;
            obj.INJ = params.INJ;
            obj.LLL = params.LLL;
            obj.DEPOP = params.DEPOP;
            
              obj.i_ul = 2;
            
            

            obj.Tdeph_ul = params.Tdeph_ul;
            %obtain the energies of the core levels for the simulation
            hbar = Constants('hbar',{'time',obj.tch})/Constants('q0');
            obj.E0 = params.E0/hbar;
            obj.loss = params.linear_loss*100/(1/obj.lch)*ones(obj.N_pts,1);
            obj.IDX = params.IDX;
            obj.TOT_RTS = zeros(1,obj.NLVLS);
            
            obj.zUL = params.zUL ; % dipole element in units of nm
            obj.zNORM = params.zNORM;
            
            % Carrier density in 1/m^3!
            obj.Ncarriers = params.Ncarriers_cm*(100^3);
            % trace and normalization constants
            obj.TRACE_RHO = ((obj.E0*1e12*obj.Ncarriers*obj.Gamma*...
                (obj.zUL*1E-9*Constants('q0'))^2)...
                /(Constants('eps0')*obj.nTHz*Constants('c')*Constants('hbar')))...
                /(1/(obj.lch*obj.tch));
            
            obj.NORM_FACTOR_DM = obj.zUL/obj.zNORM;
            obj.NORM_FACTOR_FIELD = obj.TRACE_RHO*obj.zNORM./obj.zUL;
            
            obj.global_idx_rest = setdiff(1:obj.NLVLS,[obj.INJ, ...
                obj.ULL, obj.LLL, obj.DEPOP,obj.IGNORELEVEL]);
            obj.N_rest = length(obj.global_idx_rest);
            
            obj.DEPH = zeros(3);
            obj.W = zeros(obj.NLVLS,obj.NLVLS);
            
            % population densities
            obj.rho_u = 1/3*ones(obj.N_pts,1);
            obj.rho_u_solver = MS(5,obj.N_pts,[],obj.rho_u);
            
            obj.rho_i = 1/3*ones(obj.N_pts,1);
            obj.rho_i_solver = MS(5,obj.N_pts,[],obj.rho_i);
            
            obj.rho_l =  1 - obj.rho_u-obj.rho_u;
            obj.rho_l_solver = MS(5,obj.N_pts,[],obj.rho_l);
            
            % polarizations
            % excited ground forward and backward
            obj.eta_ul = 1e-10*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.eta_ul_t = 0*obj.eta_ul;
            obj.eta_ul_solver = MS(5,obj.N_pts,[],obj.eta_ul);
            
            
            obj.rRES = cell(1,obj.N_rest);
            obj.rRES_t = cell(1,obj.N_rest);
            obj.rRES_solver = cell(1,obj.N_rest);
            
            
            for i = 1:obj.N_rest
                obj.rRES{i} = zeros(obj.N_pts,1);
                obj.rRES_t{i} = zeros(obj.N_pts,1);
                obj.rRES_solver{i} = MS(5,obj.N_pts,[],obj.rRES{i});
            end
            
            %finally calculate the scattering matrix elements and the
            %energies.
            obj.prep_data(params.Wmtx,params.HTB)
        end
        
        
        
        function [] = prep_data(obj,W,HTB_ev)
           
            INJ = obj.INJ; ULL = obj.ULL; LLL = obj.LLL; DEPOP = obj.DEPOP;
            % just copy the levels!
            for i_ = 1:obj.NLVLS
                for j_ = 1:obj.NLVLS
                    if i_ == j_
                        continue
                    end
                    obj.W(i_,j_) = W(j_+(i_-1)*obj.NLVLS);
                end
            end
            
            % ELIMINATE DEPOP level (periodic boundary conditions)
            for i_ = setdiff(1:obj.NLVLS, [INJ,DEPOP])
                obj.W(i_,INJ) = obj.W(i_,INJ) + obj.W(i_,DEPOP);
                obj.W(i_,DEPOP) = 0.;
            end
            % no backscattering from the DEPOP level;
            obj.W(DEPOP,:) = 0.;
            
            
            
            obj.TOT_RTS(INJ) = 0*obj.TOT_RTS(INJ);
            for el = setdiff(1:obj.NLVLS,[INJ,DEPOP])
                obj.TOT_RTS(INJ) = obj.TOT_RTS(INJ)+obj.W(INJ,el);
            end
            
            for lvl = setdiff(1:obj.NLVLS,[INJ,DEPOP,obj.IGNORELEVEL])
                obj.TOT_RTS(lvl) = 0* obj.TOT_RTS(lvl);
                for el = 1:obj.NLVLS
                    if el ~=lvl
                        obj.TOT_RTS(lvl) = obj.TOT_RTS(lvl) +  obj.W(lvl,el);
                    end
                end
            end
            
            % % % % Added Pure Dephasing
        
            obj.DEPH(obj.i_ul) = 1/2*(obj.TOT_RTS(obj.ULL)+...
                obj.TOT_RTS(obj.LLL))+1/obj.Tdeph_ul;
         
            
            %obtain the energies of the core levels for the simulation
            hbar = Constants('hbar',{'time',obj.tch})/Constants('q0');
            
            E_U = HTB_ev(ULL+(ULL-1)*obj.NLVLS)/hbar;
            E_L = HTB_ev(LLL+(LLL-1)*obj.NLVLS)/hbar; % in rad/ps
            
            %%%% TUNNELING transition is INJ->ULL
            
            obj.E_UL = E_U-E_L; %rad/ps  3->2 transition freq (optical transition)
            obj.dE_UL = +1i*(obj.E0 - obj.E_UL) - obj.DEPH(obj.i_ul); %
            obj.get_steady_state;
            
        end
        
        function [] = propagate(obj,FIELD,dt)
            BACKSCATTERING = true;
            FIELD = obj.NORM_FACTOR_DM.*FIELD;
            W = obj.W;
            INJ = obj.INJ;
            ULL = obj.ULL;
            LLL = obj.LLL;
            DEPOP = obj.DEPOP;
            
            %time resolved dipole moment ratio.
            %%%% POPULATIONS
            obj.rho_i_t = W(ULL,INJ).*obj.rho_u + ...
                W(LLL,INJ).*obj.rho_l - obj.TOT_RTS(INJ).*obj.rho_i;
            if(BACKSCATTERING)
                obj.rho_i_t = obj.rho_i_t-(W(DEPOP,ULL)+W(DEPOP,LLL))*obj.rho_i;
            end
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.rho_i_t = obj.rho_i_t + W(p_glob_idx,INJ).*obj.rRES{p};
                if(BACKSCATTERING)
                    obj.rho_i_t = obj.rho_i_t-W(DEPOP,p_glob_idx)*obj.rho_i;
                end
            end
            obj.rho_i_solver.make_step(obj.rho_i_t,dt);
            
            lmInteraction = conj(FIELD).*obj.eta_ul;
            obj.rho_u_t = 1i/2*(lmInteraction-conj(lmInteraction)) + ...
                obj.rho_i.*W(INJ,ULL) + obj.rho_l.*W(LLL,ULL) - obj.TOT_RTS(ULL).*obj.rho_u;
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.rho_u_t = obj.rho_u_t + W(p_glob_idx,ULL).*obj.rRES{p};
            end
            if(BACKSCATTERING)
                obj.rho_u_t = obj.rho_u_t + W(DEPOP,ULL)*obj.rho_i;
            end
            obj.rho_u_solver.make_step(obj.rho_u_t,dt);
            
            obj.rho_l_t = -1i/2*(lmInteraction-conj(lmInteraction)) + obj.rho_i.*W(INJ,LLL)+ obj.rho_u.*W(ULL,LLL) -...
                obj.TOT_RTS(LLL).*obj.rho_l;
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.rho_l_t = obj.rho_l_t + W(p_glob_idx,LLL).*obj.rRES{p};
            end
            if(BACKSCATTERING)
                obj.rho_l_t = obj.rho_l_t + W(DEPOP,LLL)*obj.rho_i;
            end
            obj.rho_l_solver.make_step(obj.rho_l_t,dt);
            
            
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.rRES_t = W(INJ,p_glob_idx).*obj.rho_i+W(ULL,p_glob_idx).*obj.rho_u+ ...
                    W(LLL,p_glob_idx).*obj.rho_l - obj.TOT_RTS(p_glob_idx).*obj.rRES{p};
                %add the rate equations part!
                for j =1:obj.N_rest % nr
                    if j~= p
                        j_glob_idx = obj.global_idx_rest(j);
                        obj.rRES_t = obj.rRES_t + W(j_glob_idx,p_glob_idx).*obj.rRES{j};
                    end
                end
                if(BACKSCATTERING)
                    obj.rRES_t = obj.rRES_t + W(DEPOP,p_glob_idx)*obj.rho_i;
                end
                obj.rRES_solver{p}.make_step(obj.rRES_t,dt);
            end
            
            
            %%%% COHERENCES
                
            obj.eta_ul_t = obj.dE_UL.*obj.eta_ul + 1i/2.*FIELD.*(obj.rho_u-obj.rho_l);
            obj.eta_ul_solver.make_step(obj.eta_ul_t,dt);
            
        
            
        end
        
        function [trace] = get_avg_trace(obj)
            trace = obj.rho_i+obj.rho_u+obj.rho_l;
            for p=1:obj.N_rest
                trace = trace + obj.rRES{p};
            end
            trace = mean(trace);
            
        end
        
        function [] = update_state(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%
            obj.rho_i = obj.rho_i_solver.get_latest_solution();
            obj.rho_u = obj.rho_u_solver.get_latest_solution();
            obj.rho_l = obj.rho_l_solver.get_latest_solution();
            for p = 1:obj.N_rest
                obj.rRES{p} = obj.rRES_solver{p}.get_latest_solution();
            end
            obj.eta_ul = obj.eta_ul_solver.get_latest_solution();
            
            
        end
        
        function [P,P_t,LOSSES] = get_polarization_and_losses(obj)
            % additionally the polarizaiton should be multiplied by -1i*c/n
            
            % forward wave polarization
            P = obj.NORM_FACTOR_FIELD.*obj.eta_ul;
            P_t = obj.NORM_FACTOR_FIELD.*obj.eta_ul_t;
            % losses
            LOSSES = obj.loss;
        end
        
        function [J_TL] = get_current_density(obj,params)
            %%% returns the section's current density in units A/mm^2,
            %%% calculated based on the periodic rate equations' formula.
            
            %current density will be calculated with SHB included!
            r11 = obj.rho_i;
            r33 = obj.rho_u;
            r22 = obj.rho_l;
            
            rates =   r11.*obj.W(obj.INJ,obj.DEPOP) +  r33.*obj.W(obj.ULL,obj.DEPOP)+...
                r22.*obj.W(obj.LLL,obj.DEPOP);
            for p=1:obj.N_rest
                rates = rates + obj.rRES{p}.*obj.W(obj.global_idx_rest(p),obj.DEPOP);
            end
            
            J_TL =  (Constants('q0')*(obj.Lp*1E-9)*obj.Ncarriers*1E12)/1E6*rates; %in A/mm^2
            
            % this is the leading expression of the tunnelign current.
            %             rho_iu = obj.rho_iu+;
            %             rates2 =   1i*obj.OMEGA.*(rho_iu-conj(rho_iu));
        end
        
       
        function steady_state = get_steady_state(obj)
            %%%%%
            
            % active levels
            active_levels = setdiff(1:obj.NLVLS,[obj.DEPOP,obj.IGNORELEVEL]);
            GAMMA = transpose(obj.W);

            GAMMA = GAMMA(active_levels,active_levels);
            
            for i=1:length(GAMMA)
                GAMMA(i,i) = 0.
                GAMMA(i,i) = -sum(GAMMA(:,i));
            end
            
            NLV = length(GAMMA);
            % level index to "eliminate" by the boundary condition
            % sum_i rho_ii = 1
            el_ = NLV;
            obj.S = zeros(NLV-1,NLV-1);
            ids = setdiff(1:NLV,el_);
            obj.b = GAMMA(ids,el_);
            for idx_i =1:length(ids)
                for idx_j = 1:length(ids)
                    obj.S(idx_i,idx_j) = GAMMA(ids(idx_i),ids(idx_j))-GAMMA(ids(idx_i),el_);
                end
            end
            
            % THIS STEADY STATE IS ALSO INFLUENCED BY THE TUNNELING
            % COUPLING!
            steady_state = -inv(obj.S)*obj.b;
            obj.steady_state = [steady_state; 1- sum(steady_state)];
            
            [P,d] = eigs(obj.S,length(obj.S));
            
            P_inv = inv(P);
            Tmtx = -inv(d);
            
        end
        
        function d_th = calc_threshold_inv(obj,params_gain,params_abs)
            
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
            
            sigma_g = (1./obj.DEPH(obj.i_ul))*(1e-12)*params_gain.Gamma*E0*1E12*mu_g^2/ ...
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
            obj.d_th = d_th;
        end
        
        
        function T_GR = calculate_gain_recovery(obj)
            
            pops_ = [obj.rho_i,obj.rho_u,obj.rho_l];
            
%             do not include this level since it is the EOM of 
%             this state that is used for the normalization condition; 

%             for pop_ctr = 1: length(obj.rRES)
%                 pops_ = [pops_, obj.rRES{pop_ctr}];
%             end
            
            n = transpose(pops_);
            n_tilde = 0*n;
            
            T_GR = 0;
            [P,d] = eigs(obj.S);
            
            P_inv = inv(P);
            for j = 1:size(n,2)
                n_tilde(:,j) = P_inv*n(:,j);
            end
            
            Tmtx = -inv(d);
            
            inv_eq = obj.steady_state(obj.ULL)-obj.steady_state(obj.LLL); 
            n_tilde_eq = Tmtx*P_inv*obj.b;
            inv_vtr = (obj.d_th - inv_eq)*ones(1,size(n,2));
            function rhs = inversion_recovery(tau)
                rhs = inv_vtr;
                for lvl =1:size(n,1) 
                    rhs(:) = rhs(:) - (P(obj.ULL,lvl)-P(obj.LLL,lvl))*(n_tilde(lvl,:) - n_tilde_eq(lvl))*exp(-tau(:)./Tmtx(lvl,lvl));
                end
            end
            func_ptr = @inversion_recovery;
            options = optimoptions('fsolve');
%             options.TolFun = 1e-3; 
%             options.MaxFunEvals = 100; % otherwise slow! 
            T_GR = fsolve(func_ptr,zeros(size(n(1,:))),options); 
            
        end
        
        
        
        
    end
end

