function dat = TL_voltage(rlgc,dat)
%%%%% TRANSMISSION LINE EQUATIONS %%%%%%%%%%%%%
% The voltage source (Vs) with wire impedance (Zin) is replace with its Northon
% equivalent circuit, which consists of a current source (Vs/Zin) and a
% resist Zin in parallel.

%Definition of constants
dat.dx = 1;        %unit: mm
dat.dt = 1;        %unit: s
dat.Zin = 50;      %unit: ohm
dat.Vs = 9;        %unit: V

rlgc.C1 = 1.5e-12;  %unit: F/mm
rlgc.L1 = 1.6e-10;  %unit: H/mm
rlgc.C2 = 1.5e-12;  %unit: F/mm
rlgc.L2 = 1.6e-10;  %unit: H/mm

dat.w1 = 50e-3;     %unit: mm
dat.h1 = 10e-3;     %unit: mm
dat.w2 = 50e-3;     %unit: mm
dat.h2 = 10e-3;     %unit: mm
%Initial conditions
%At t=0, the voltage on lines were set to be equal with half of the source,
%and no current flows to lines. Due to tunneling of electrons (J) which is 
%relevant to position, the voltage as well as current along the lines 
%will be changed until equilibrium.
dat.v_1TL(1:end) = dat.Vs/2;    % V(0)    V(1)....   V(M)
dat.i_1TL(1:end) = 0;           % I(1/2)  I(3/2).... I(M+1/2)
dat.J_1TL(1:end) = 0;           % J(0)    J(1).....  J(M)      A/mm^2

dat.v_2TL(1:end) = dat.Vs/2;    % V(0)    V(1)....   V(M)
dat.i_2TL(1:end) = 0;           % I(1/2)  I(3/2).... I(M+1/2)
dat.J_2TL(1:end) = 0;           % J(0)    J(1).....  J(M)      A/mm^2

%Boundary conditions
%right side open circuit,image current
dat.i_1TL(end) = -dat.i_1TL(end-1);
dat.v_1TL(end) = dat.v_1TL(end)-dat.dt/rlgc.C1/dat.dx*(-2*dat.i_1TL(end-1)-...
                    -dat.w1*dat.dx*dat.J_1TL(end));

dat.i_2TL(end) = -dat.i_2TL(end-1);
dat.v_2TL(end) = dat.v_2TL(end)-dat.dt/rlgc.C/dat.dx*(-2*dat.i_2TL(end-1)-...
                    -dat.w2*dat.dx*dat.J_2TL(end));


%Voltage and current updating
        % Section 1 §
temp_v1 = dat.v_1TL(1); % store the old voltage on section 1
dat.v_1TL(1) = 1/(rlgc.C1*dat.dx/dat.dt/2+1/2/dat.Zin)*...
                    ((rlgc.C1*dat.dx/dat.dt/2-1/2/dat.Zin)*...
                    dat.v_1TL(1)-dat.i_1TL(1)-dat.width*dat.dx/2*...
                    dat.J_1TL(1)+(dat.Vs-dat.v_2TL(1))/dat.Zin);

dat.v_1TL(2:end-1) = dat.v_1TL(2:end-1)-dat.dt/rlgc.C1/dat.dx*(dat.i_1TL(2:end-1)-...
                    dat.i_1TL(1:end-2)-dat.w1*dat.dx*dat.J_1TL(2:end-1));

dat.i_1TL(1:end-1) = dat.i_1TL(1:end-1)-dat.dt/rlgc.L1/dat.dx*(dat.v_1TL(2:end)...
                    -dat.v_1TL(1:end-1));

        % Section 2 §
dat.v_2TL(1) = dat.Vs-dat.v_1TL(1)-dat.Zin*(dat.i_1TL(1)+dat.w1*dat.dx*dat.J_1TL(1)+...
                rlgc.C1*dat.dx*dat.dt/2*(dat.v_1TL(1)-temp_v1));

dat.v_2TL(2:end-1) = dat.v_2TL(2:end-1)-dat.dt/rlgc.C2/dat.dx*(dat.i_2TL(2:end-1)-...
                    dat.i_2TL(1:end-2)-dat.w2*dat.dx*dat.J_2TL(2:end-1));

dat.i_2TL(1:end-1) = dat.i_2TL(1:end-1)-dat.dt/rlgc.L2/dat.dx*(dat.v_2TL(2:end)...
                    -dat.v_2TL(1:end-1));

                 
                
end