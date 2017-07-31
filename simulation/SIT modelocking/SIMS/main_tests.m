Ls = 1.1:.1:1.5; % length of the gain material 
Ns = linspace(4e16,5.5e16,5);

% vary the ratio of saturation 
p = 1.2; % 1.2 times above thold; 
T2_abs = .5:.25:1.5; 

for T2_select = T2_abs; 

    params_abs = input_parser('SIMFILES/ABSORBER_INPUTS.sim');
    params_gain = input_parser('SIMFILES/GAIN_INPUTS.sim');
    params_abs.T2 = T2_select;
    model = '2LVL';
    

    
    % modify params_gain if needed?
    prefix = 'PML_2LVL_FP'
    T1g = params_gain.T_1; 
    T2g = params_gain.T_2;
    T1gstr = strsplit(num2str(T1g,'%.2f'),'.'); 
    T1gstr = [T1gstr{1} 'p' T1gstr{2}];
    T2gstr = strsplit(num2str(T2g,'%.2f'),'.'); 
    T2gstr = [T2gstr{1} 'p' T2gstr{2}];
    pstr = strsplit(num2str(p,'%.2f'),'.'); pstr = [pstr{1} 'p' pstr{2}];

    
    T1a = params_abs.T_1; T2a = params_abs.T_2;

    dipole_r = params_abs.zUL/params_gain.zUL;
    
    T1astr = strsplit(num2str(T1a,'%.2f'),'.');
    T1astr = [T1astr{1} 'p' T1astr{2}];
    T2astr = strsplit(num2str(T2a,'%.2f'),'.'); T2astr = [T2astr{1} 'p' T2astr{2}];
    

    dipolestr = strsplit(num2str(dipole_r,'%.2f'),'.'); dipolestr = [dipolestr{1} 'p' dipolestr{2}];
    
    savename = ['SIMRES/' prefix ';T1a=' T1astr ',T2a=' T2astr ';T1g=' T1gstr ',T2g=' T2gstr ';' ...
         pstr ';dipole_ratio=' dipolestr];
%     passive_modelocking_FP(params_gain,params_abs,p,savename,model);
     passive_modelocking(params_gain,params_abs,p,savename,model);
end