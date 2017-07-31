function [ const ] = Constants(name,varargin)
% This function returns a struct or constants or a single constant from 
% a predefined list with physical constants for our simulations. 

% Input: 
% the user gives a string , naming the constant, or an empty string if 
% he/she wants a full list of constants. Available constants till now are: 
%
%           Description           Name             Base Units: 
%    speed of light in vacuum:    'c'                m/s
%    Electric Constant            'eps0'             F/m
%    Planck constnat/2*pi         'hbar'             J.s
%    permeability                 'mu0'              H/m 
%    Electron mass                'me'               kg
%    Boltzmann constant           'kb'               J/K
%    elementary charge            'q0'               Coulombs

% The user has also the option to recover any of the above mentioned
% constants in a custom defined system of units, by simply giving the
% characteristic unit for the new system. For example, if he/she wishes to
% get Plancks constant in J-ps (i.e. characteristic time is picosecond) he
% has to call Constants('hbar',{'time',1e-12}); where 'time' is the quantity
% the user wishes to change the units of and 1e-12 is the conversion
% number that converts quantities from the old system to the new system:
% i.e. 1ps = 1e-12 seconds. User can change the system of units 
% by specifying characteristic values only for one or more of he fundamental 
% quantities (time, length, mass , and temperature) and not for any
% derived quantity. Notice that ONLY explicit conversions are supported. 
% for example if one calls Constants('hbar',{'time',1e-12}) the function
% will return hbar in units Joules-ps. instead of kg m^2/ps ! 
% 


% constants in SI units. 
c = 2.99792458e+8;        % Speed of light in vacuum, m/s
eps0 = 8.854187817e-12;   % Electric constant, F/m (or Coulomb^2/Joules) 
mu0  = 4*pi*1E-7;         % permeability in H/m 
hbar = 1.05457266e-34;    % Planck constant/(2*pi), J-s 
q0 = 1.60217733e-19;     % Elementary charge Coulombs
me = 9.1093897e-31;       % Electron mass, kg
kb = 1.380658e-23;         % Boltzmann constant, J/Kelvin


 for i = 1:length(varargin)
       argi = varargin{i};
       
       if(strcmp(argi{1},'time'))
        tch =argi{2};
        c = c/(1/tch); 
        hbar = hbar*1/tch; %Joule - tch
       end       
       if(strcmp(argi{1},'length'))
        lch = argi{2};
        c = c/lch;
        eps0 = eps0/1/lch; % F/lch
        mu0 = mu0/1/lch; %H/lch
       end
       
       if(strcmp(argi{1},'mass'))
        mch = arg{2};
        me = me*1/mch; 
       end
       if(strcmp(argi{1},'temperature')) %only boltzmann constant has temperature in it
        tempch = argi{2}; 
        kb = kb/(1/tempch);
       end
 
 end

if(strcmp(name,''))
 const.c =c ;        % Speed of light in vacuum, m/s
 const.eps0 = eps0;   % Electric constant, F/m (or Coulomb^2/Joules) 
 const.mu0 = mu0;     % permeability in H/m 
 const.hbar = hbar;    % Planck constant/(2*pi), J-s 
 const.q0 = q0;     % Elementary charge Coulombs
 const.me = me;       % Electron mass, kg
 const.kb =kb;  
   
else
    if(strcmp(name,'c'))
        const = c; 
    end
    if(strcmp(name,'eps0'))
        const = eps0;
    end
    if(strcmp(name,'mu0'))
        const = mu0;
    end
    if(strcmp(name,'hbar'))
        const = hbar; 
    end
    if(strcmp(name,'me'))
        const = me; 
    end
    if(strcmp(name,'kb'))
        const = kb; 
    end
    if(strcmp(name,'q0'))
        const = q0; 
    end
      
end
    
    
end

