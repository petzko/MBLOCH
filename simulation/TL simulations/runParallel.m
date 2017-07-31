clear
foldername = 'sim_scenarios'
interpDataFile = 'fitted_data_OPTICA.mat';
fold_cont = dir(foldername); 
pattern = 'sec\w*'
simfilenames = {};
savefilenames = {};
ctr = 1; 

for f_idx = 1:length(fold_cont)
   f = fold_cont(f_idx);
   fname =  f.name;
   startidx = regexp(fname,pattern,'match');
   if length(startidx) >0 
       fullname = [foldername,'/',fname]
       
       simset_ = input_parser(fullname);
       modA= simset_.modA; 
       strmodA = strsplit(num2str(modA,'%.2f'),'.');
       strmodA = [strmodA{1},'p',strmodA{2}];
       
       modF = simset_.modF*1e3; 
       strmodF = strsplit(num2str(modF,'%.2f'),'.');
       strmodF = [strmodF{1},'p',strmodF{2}];
       
       bias = simset_.bias;
       strbias = strsplit(num2str(bias,'%.2f'),'.');
       strbias = [strbias{1} 'p'   strbias{2}];
       
       savename = ['ACmodelocking;modAmpl=' strmodA ';modFreq=' strmodF 'GHz;init_bias=' strbias 'kVpcm;'];
       
       savefilenames{ctr} = ['sim_results/' savename];
       simfilenames{ctr} = fullname
       ctr = ctr+1;
   
   end
end


Ni = length(simfilenames); 
workers = Ni;

spmd(workers)
    tidx = labindex;    
     acmodelocking(simfilenames{tidx},interpDataFile, savefilenames{tidx},false);
end
