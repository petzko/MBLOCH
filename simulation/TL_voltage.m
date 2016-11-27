function dat = TL_voltage(rlgc,dat)
%%%%% TRANSMISSION LINE EQUATIONS %%%%%%%%%%%%%
% The voltage source (Vs) with wire impedance (Zin) is replace with its Northon
% equivalent circuit, which consists of a current source (Vs/Zin) and a
% resist Zin in parallel.

%Definition of constants
rlgc.C = 1.5;  %unit: pF/mm
rlgc.L = 0.16; %unit: nH/mm
dat.Zin = 50;  %unit: ohm
dat.Vs = 9;    %unit: V

%Initial conditions
%At t=0, the voltage on the line was set to be equal with the source, and
%no current flows to the line. Due to tunneling of electrons (J) which is 
%relevant to the position, the voltage as well as current along the line 
%will be changed until to the equilibrium.
dat.v_TL(1:end) = dat.Vs; % V(0)    V(1)....   V(M)
dat.i_TL(1:end) = 0;      % I(1/2)  I(3/2).... I(M+1/2)

%Boundary conditions
dat.i_TL(end) = 0; % right side open circuit

%Voltage and current updating
dat.v_TL(1:end) = 1/(rlgc.C*dat.dx/dat.dt+1/2/dat.Zin)*...
                    ((rlgc.C*dat.dx/dat.dt-1/2/dat.Zin)*...
                    dat.v_TL(1:end)-dat.i_TL(1:end)+dat.J(1:end)+dat.Vs/dat.Zin);

dat.i_TL(1:end-1) = dat.i_TL(1:end-1)-dat.dt/rlgc.L/dat.dx*(dat.v_TL(2:end)...
                    -dat.v_TL(1:end-1));
                