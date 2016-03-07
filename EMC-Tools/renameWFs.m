wfreadOFFSET = 50; 
WF_grid = dlmread('input1.txt','\t',wfreadOFFSET,0);
fid = fopen('PSI.dat','w');
fprintf(fid, '       %i  %.15E \r\n',WF_grid');
fclose(fid);