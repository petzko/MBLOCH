clear


sim_scanarios_foldername = 'sim_scenarios' ;
save_files_folder_name = 'IVCHAR' ; 

interpDataFile = 'fitted_data_OPTICA.mat';
fold_cont = dir(sim_scanarios_foldername); 
simfilename = [sim_scanarios_foldername '/CUSTOMSET.set'];
bias_vec = [6.0:0.2:14.];
sim_settings_array = {}; 
savefilenames = {};

ctr = 1; 
for  idx = 1:length(bias_vec)
       simset_ = input_parser(simfilename);
       simset_.bias = bias_vec(idx);

       strbias = strsplit(num2str(bias_vec(idx),'%.2f'),'.');
       strbias = [strbias{1} 'p'   strbias{2}];
       
       savename = ['IVLchar;at-bias=' strbias 'kVpcm;'];
       
       savefilenames{idx} = [save_files_folder_name '/' savename];
       sim_settings_array{idx} = simset_
        
end


workers = min(16,round(length(bias_vec)/4)); 
% split the job
dN = ceil(length(bias_vec)/workers);
jobs_ =cell(workers,1);
for tid =0:(workers-1)
    jobs_{tid+1}.a = tid*dN+1;
    jobs_{tid+1}.b = min((tid+1)*dN,length(bias_vec));
end

spmd(workers)
    tidx = labindex;    
    for idx =jobs_{tidx}.a:jobs_{tidx}.b
        acmodelocking('DUMMY',interpDataFile, [savefilenames{idx},'byworker=' num2str(tidx)],false,'simsettings',sim_settings_array{idx},'simrt',100);
    end
end

