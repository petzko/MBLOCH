classdef DM_MODEL_PSEUDO_3_LVL_RWA_FP < handle
    properties
        name;
        
        % excited state
        r330;  r330_solver;   r330_t;
        r33p;  r33p_solver;   r33p_t;
        
        % lower laser state
        r220; r220_solver;  r220_t;
        r22p; r22p_solver;  r22p_t;
        
        % depopulation state (RT-coupled  to the excited state)
        r110;  r110_solver;  r110_t;
        r11p;  r11p_solver;  r11p_t;
        
        
        
        % rest of the levels~
        rRES; rRES_t; rRES_solver;
        rRES_SHB; rRES_SHB_t; rRES_SHB_solver;
        % coherences backward and forward components.
        n32p; n32m;
        
        n32p_t;   n32m_t;
        n32p_solver;  n32m_solver;
        
        
        % Energies in units (1/time)
        E12, E32, E13, O13;
        dE12,dE32,dE13;
        % central frequency
        E0;hbar;
        % period length
        Lp;
        % spatial hole burning
        shb;
        
        % Indices of the Excited, Spin and the Ground state as well as the
        % depopulation level
        ULL, INJ, LLL, DEPOP,IGNORELEVEL;
        
        % total number of levels
        NLVLS;
        global_idx_rest;
        N_rest;
        
       
        % scattering rates smatrix
        W;
        % 1/lifetime array!
        TOT_RTS;
        DEPH;
        Tdeph_32;
        IDX;
        
        % current density
        J_0; J_2;
        
        
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
        Overlap,zUL,zNORM;
        Ncarriers;
        loss,nTHz;
        tch,lch,diffusion; % characteristic time and length
    end
    
    
    methods
        
        function obj  = DM_MODEL_PSEUDO_3_LVL_RWA_FP(params)
            
            obj.name = params.name;
            %characteristic length and time
            obj.tch = params.tch;
            obj.lch = params.lch;
            % ALLlevels
            obj.NLVLS = params.NLVLS;
            obj.Lp = params.Lp;
            % nomalization constant params
            obj.Overlap = params.Overlap;  % overlap factor
            obj.zUL = params.zUL ; % dipole element in units of nm
            obj.nTHz = params.nTHz;
            obj.N_pts = params.N_pts;
            
            obj.ULL = params.ULL;
            obj.INJ = params.INJ;
            obj.LLL = params.LLL;
            obj.DEPOP = params.DEPOP;
            obj.IGNORELEVEL = params.IGNORELEVEL;
            
            obj.zNORM = params.zNORM;
            obj.Tdeph_32 = params.Tdeph_32;
            
            obj.E0 = params.E0;
            obj.hbar = Constants('hbar',{'time',obj.tch})/Constants('q0');

            obj.loss = params.linear_loss*100/(1/obj.lch)*ones(obj.N_pts,1);
            obj.IDX = params.IDX;
            obj.diffusion = 0*4*obj.E0/params.c^2*params.D*10^2/(1/params.tch);
            obj.TOT_RTS = zeros(obj.N_pts,obj.NLVLS);
            
            obj.shb = params.shb;
            
            % Carrier density in 1/m^3!
            obj.Ncarriers = params.Ncarriers_cm*(100^3);
            % trace and normalization constants
            obj.TRACE_RHO = ((obj.E0*1e12*obj.Ncarriers*obj.Overlap*...
                (obj.zUL*1E-9*Constants('q0'))^2)...
                /(Constants('eps0')*obj.nTHz*Constants('c')*Constants('hbar')))...
                /(1/(obj.lch*obj.tch));
            
            
            obj.global_idx_rest = setdiff(1:obj.NLVLS,[obj.INJ, ...
                obj.ULL, obj.LLL, obj.DEPOP,obj.IGNORELEVEL]);
            obj.N_rest = length(obj.global_idx_rest);
            
            obj.DEPH = zeros(obj.N_pts,1);
            obj.W = zeros(obj.N_pts,obj.NLVLS,obj.NLVLS);
            
            % population densities
            obj.r330 = 0*ones(obj.N_pts,1);
            obj.r330_solver = MS(5,obj.N_pts,[],obj.r330);
            
            obj.r110 = 1/3*ones(obj.N_pts,1);
            obj.r110_solver = MS(5,obj.N_pts,[],obj.r110);
            
            obj.r220 =  1 - obj.r330-obj.r110;
            obj.r220_solver = MS(5,obj.N_pts,[],obj.r220);
            
            
            % shb terms
            obj.r33p = zeros(obj.N_pts,1);
            obj.r33p_solver = MS(5,obj.N_pts,[],obj.r33p);
            obj.r22p =  zeros(obj.N_pts,1);
            obj.r22p_solver = MS(5,obj.N_pts,[],obj.r22p);
            obj.r11p = zeros(obj.N_pts,1);
            obj.r11p_solver = MS(5,obj.N_pts,[],obj.r11p);
            
            % polarizations
            % excited ground forward and backward
            obj.n32p = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.n32p_t = 0*obj.n32p;
            obj.n32p_solver = MS(5,obj.N_pts,[],obj.n32p);
            
            obj.n32m = 1e-15*((rand(obj.N_pts,1)-1/2)+ ...
                1i*(rand(obj.N_pts,1)-1/2));
            obj.n32m_t = 0*obj.n32m;
            obj.n32m_solver = MS(5,obj.N_pts,[],obj.n32m);
            
            obj.rRES = cell(1,obj.N_rest); obj.rRES_SHB = cell(1,obj.N_rest);
            obj.rRES_t = cell(1,obj.N_rest); obj.rRES_SHB_t = cell(1,obj.N_rest);
            obj.rRES_solver = cell(1,obj.N_rest);
            obj.rRES_SHB_solver = cell(1,obj.N_rest);
            
            for i = 1:obj.N_rest
                obj.rRES{i} = zeros(obj.N_pts,1);
                obj.rRES_t{i} = zeros(obj.N_pts,1);
                obj.rRES_solver{i} = MS(5,obj.N_pts,[],obj.rRES{i});
                obj.rRES_SHB{i} = zeros(obj.N_pts,1);
                obj.rRES_SHB_t{i} = zeros(obj.N_pts,1);
                obj.rRES_SHB_solver{i} = MS(5,obj.N_pts,[], obj.rRES_SHB{i});
            end
            
            obj.J_0 = zeros(obj.N_pts,1);
            obj.J_2 = zeros(obj.N_pts,1);
            
        end
        
        function [] = propagate(obj,U,V,dt)
            
            U = obj.NORM_FACTOR_DM.*U(obj.IDX);
            V = obj.NORM_FACTOR_DM.*V(obj.IDX);
            
            W = obj.W; INJ = obj.INJ; ULL = obj.ULL; LLL = obj.LLL;  DEPOP = obj.DEPOP;
            
            %time resolved dipole moment ratio.
            %%%% POPULATIONS
            
            Temp = 50; % Kelvin
            beta_ps =(Constants('kb')*Temp); %J
            beta_ps = (beta_ps/Constants('q0')/obj.hbar)^(-1); %ps/rad
            TEMPDEP = exp(-beta_ps.*abs(obj.E13));
            gamma2 = 2*obj.DEPH;
            
            PRECONST = obj.O13.^2.*gamma2./(obj.E13.^2+gamma2.^2); % todo implement
            obj.J_0 = PRECONST.*(heaviside(obj.E13).*(obj.r110-TEMPDEP.*obj.r330)+ ...
                heaviside(-obj.E13).*(obj.r110.*TEMPDEP - obj.r330));
            
            lmInteraction = (conj(U).*obj.n32p + conj(V).*obj.n32m);
            obj.r330_t = obj.J_0 + 1i/2*(lmInteraction-conj(lmInteraction)) + ...
                obj.r110.*W(:,DEPOP,ULL) + obj.r220.*W(:,LLL,ULL) - obj.TOT_RTS(:,ULL).*obj.r330;
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.r330_t = obj.r330_t + W(:,p_glob_idx,ULL).*obj.rRES{p};
            end
            obj.r330_solver.make_step(obj.r330_t,dt);
            
            obj.r220_t = -1i/2*(lmInteraction-conj(lmInteraction)) + ...
                obj.r110.*W(:,DEPOP,LLL)+ obj.r330.*W(:,ULL,LLL) -...
                obj.TOT_RTS(:,LLL).*obj.r220;
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.r220_t = obj.r220_t + W(:,p_glob_idx,LLL).*obj.rRES{p};
            end
            obj.r220_solver.make_step(obj.r220_t,dt);
            
            obj.r110_t = - obj.J_0 +  W(:,ULL,DEPOP).*obj.r330 + ...
                W(:,LLL,DEPOP).*obj.r220 - obj.TOT_RTS(:,DEPOP).*obj.r110;
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.r110_t = obj.r110_t + W(:,p_glob_idx,DEPOP).*obj.rRES{p};
            end
            obj.r110_solver.make_step(obj.r110_t,dt);
            for p = 1:obj.N_rest
                p_glob_idx = obj.global_idx_rest(p);
                obj.rRES_t = W(:,ULL,p_glob_idx).*obj.r330 + W(:,LLL,p_glob_idx).*obj.r220 + ...
                   W(:,DEPOP,p_glob_idx).*obj.r110 - obj.TOT_RTS(:,p_glob_idx).*obj.rRES{p};
                %add the rate equations part!
                for j =1:obj.N_rest % nr
                    if j~= p
                        j_glob_idx = obj.global_idx_rest(j);
                        obj.rRES_t = obj.rRES_t + W(:,j_glob_idx,p_glob_idx).*obj.rRES{j};
                    end
                end
                obj.rRES_solver{p}.make_step(obj.rRES_t,dt);
            end
            
            
            %%%% COHERENCES
            
            obj.n32p_t = obj.dE32.*obj.n32p + 1i/2.*(U.*(obj.r330-obj.r220) + V.*(obj.r33p-obj.r22p));
            obj.n32p_solver.make_step(obj.n32p_t,dt);
            
            obj.n32m_t = obj.dE32.*obj.n32m + 1i/2.*(V.*(obj.r330-obj.r220) + U.*conj(obj.r33p-obj.r22p));
            obj.n32m_solver.make_step(obj.n32m_t,dt);
            
            
            if(obj.shb > 0 )
%                 SHB_t = 0*obj.r110;
                
                obj.J_2 = PRECONST.*(heaviside(obj.E13).*(obj.r11p-TEMPDEP.*obj.r33p)+ ...
                    heaviside(-obj.E13).*(obj.r11p.*TEMPDEP - obj.r33p));
             
                %%% r33+
                lminteraction2 = conj(V).*(obj.n32p) -(U).*conj(obj.n32m);
                obj.r33p_t = obj.J_2 + 1i/2.*lminteraction2 + ...
                     W(:,LLL,ULL).*obj.r22p + W(:,DEPOP,ULL).*obj.r11p  - (obj.TOT_RTS(:,ULL)+obj.diffusion).*obj.r33p;
                for p = 1:obj.N_rest
                    p_glob_idx = obj.global_idx_rest(p);
                    obj.r33p_t = obj.r33p_t + W(:,p_glob_idx,ULL).*obj.rRES_SHB{p};
                end
                obj.r33p_solver.make_step(obj.r33p_t,dt);
%                 SHB_t = SHB_t+obj.r33p_t;
                
                %%% r22+
                obj.r22p_t = -1i/2.*lminteraction2 + W(:,ULL,LLL).*obj.r33p + ...
                    W(:,DEPOP,LLL).*obj.r11p - (obj.TOT_RTS(:,LLL)+obj.diffusion).*obj.r22p;
                for p = 1:obj.N_rest
                    p_glob_idx = obj.global_idx_rest(p);
                    obj.r22p_t = obj.r22p_t + W(:,p_glob_idx,LLL).*obj.rRES_SHB{p};
                end
                obj.r22p_solver.make_step(obj.r22p_t,dt);
%                 SHB_t = SHB_t+obj.r22p_t;

                %%% r11+
                obj.r11p_t = -obj.J_2 + W(:,ULL,DEPOP).*obj.r33p+ ...
                    W(:,LLL,DEPOP).*obj.r22p - (obj.TOT_RTS(:,DEPOP)+obj.diffusion).*obj.r11p;
                for p = 1:obj.N_rest
                    p_glob_idx = obj.global_idx_rest(p);
                    obj.r11p_t = obj.r11p_t + W(:,p_glob_idx,DEPOP).*obj.rRES_SHB{p};
                end
                obj.r11p_solver.make_step(obj.r11p_t,dt);
%                 SHB_t = SHB_t+obj.r11p_t;

                
                %%% SHB in the rest levels
                for p = 1:obj.N_rest
                    p_glob_idx = obj.global_idx_rest(p);
                    obj.rRES_SHB_t = W(:,DEPOP,p_glob_idx).*obj.r11p+W(:,ULL,p_glob_idx).*obj.r33p+...
                        W(:,LLL,p_glob_idx).*obj.r22p - (obj.TOT_RTS(:,p_glob_idx)+obj.diffusion).*obj.rRES_SHB{p};
                    %add the rate equations part!
                    for j =1:obj.N_rest % nr
                        if j~= p
                            j_glob_idx = obj.global_idx_rest(j);
                            obj.rRES_SHB_t = obj.rRES_SHB_t + W(:,j_glob_idx,p_glob_idx).*obj.rRES_SHB{j};
                        end
                    end
                    obj.rRES_SHB_solver{p}.make_step(obj.rRES_SHB_t,dt);
%                     SHB_t = SHB_t+obj.rRES_SHB_t;
                end
                
            end
%             display(['sum shb_t' num2str(sum(SHB_t))])
            
        end
        
        function [trace_0,trace_2] = get_avg_trace(obj)
            trace_0 = obj.r110+obj.r330+obj.r220;
            trace_2 = obj.r11p+obj.r33p+obj.r22p;
            
            for p=1:obj.N_rest
                trace_0 = trace_0 + obj.rRES{p};
                trace_2 = trace_2 + obj.rRES_SHB{p};
            end
            trace_0 = mean(trace_0);
            trace_2 = abs(mean(trace_2));
            
        end
        
        function [] = interpolate(obj,v_TL,W_fit,E_fit,AC_fit,zUL_fit)
            
            %%% interpolate the structure based on the bias (i.e. voltage
            %%% per unit length) at each point of the section grid! The
            %%% bias should be normalized in units kV/mm!
            
            INJ = obj.INJ; ULL = obj.ULL; LLL = obj.LLL; DEPOP = obj.DEPOP;
            
            for i_ = 1:obj.NLVLS
                for j_ = 1:obj.NLVLS
                    if i_ == j_
                        continue
                    end
                    obj.W(:,i_,j_) = W_fit{i_,j_}(v_TL);
                end
            end
            
            %%% bi_32in i_13tting up TL params!
            obj.TOT_RTS(:,DEPOP) = 0*obj.TOT_RTS(:,DEPOP);
            for el = setdiff(1:obj.NLVLS,[INJ,DEPOP,obj.IGNORELEVEL])
                obj.TOT_RTS(:,DEPOP) = obj.TOT_RTS(:,DEPOP)+obj.W(:,DEPOP,el);
            end
            
            for lvl = setdiff(1:obj.NLVLS,[INJ,obj.IGNORELEVEL])
                obj.TOT_RTS(:,lvl) = 0* obj.TOT_RTS(:,lvl);
                for el = setdiff(1:obj.NLVLS,[lvl,INJ,obj.IGNORELEVEL])
                    obj.TOT_RTS(:,lvl) = obj.TOT_RTS(:,lvl) +  obj.W(:,lvl,el);
                end
            end
            
            % % % % Added Pure Dephasing
            
            obj.DEPH = 1/2*(obj.TOT_RTS(:,obj.ULL)+...
                obj.TOT_RTS(:,obj.LLL))+1/obj.Tdeph_32;
            
            %%%%% correctly i_13tup the energies, the anticrossings and the scattering rates...
            %obtain the energies of the core levels for the simulation
            hbar = Constants('hbar',{'time',obj.tch})/Constants('q0');
            
            E1 = E_fit{INJ}(v_TL)/hbar;
            E3 = E_fit{ULL}(v_TL)/hbar;
            E2 = E_fit{LLL}(v_TL)/hbar; % in rad/ps
            
            %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition
            %%%% is 1->3
            %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
            obj.O13 = AC_fit{1}(v_TL)/hbar; % in rad/ps
            
            obj.E13 = E1-E3; %rad/ps; 1->3 traisition freq
            obj.E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
            
            obj.dE32 = +1i*(obj.E0 - obj.E32) - obj.DEPH;
            
            %the varying dipole ratio.
            obj.zUL = zUL_fit(v_TL);
            obj.NORM_FACTOR_DM = obj.zUL/obj.zNORM;
            obj.NORM_FACTOR_FIELD = obj.TRACE_RHO*obj.zNORM./obj.zUL;
            
        end
        
        function [] = update_state(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%
            obj.r330 = obj.r330_solver.get_latest_solution();
            obj.r220 = obj.r220_solver.get_latest_solution();
            obj.r110 = obj.r110_solver.get_latest_solution();
            
            for p = 1:obj.N_rest
                obj.rRES{p} = obj.rRES_solver{p}.get_latest_solution();
            end
            
            obj.n32p = obj.n32p_solver.get_latest_solution();
            obj.n32m = obj.n32m_solver.get_latest_solution();
            
            
            if(obj.shb > 0)
                obj.r11p = obj.r11p_solver.get_latest_solution();
                obj.r33p = obj.r33p_solver.get_latest_solution();
                obj.r22p = obj.r22p_solver.get_latest_solution();
                for p = 1:obj.N_rest
                    obj.rRES_SHB{p} = obj.rRES_SHB_solver{p}.get_latest_solution();
                end
            end
        end
        
        function [P,P_t,M,M_t,LOSSES] = get_polarization_and_losses(obj)
            % additionally the polarizaiton should be multiplied by -1i*c/n
            
            % forward wave polarization
            P = obj.NORM_FACTOR_FIELD.*obj.n32p;
            P_t = obj.NORM_FACTOR_FIELD.*obj.n32p_t;
            
            % backward wave polarization
            M = obj.NORM_FACTOR_FIELD.*obj.n32m;
            M_t = obj.NORM_FACTOR_FIELD.*obj.n32m_t;
            
            % losses
            LOSSES = obj.loss;
        end
        
        function [current] = get_current_density(obj,params)
            %%% returns the section's current density in units A/mm^2,
            %%% calculated based on the periodic rate equations' formula.
            
            current = obj.J_0;
            if(obj.shb > 0)
                J_2 = obj.J_2.*exp(2*1i*params.x*obj.E0/params.c);
                current = current + J_2 + conj(J_2);
            end
            
            %current density will be calculated with SHB included!
            current =  (Constants('q0')*(obj.Lp*1E-9)*obj.Ncarriers*1E12)/1E6*current; %in A/mm^2
        end
    end
end

