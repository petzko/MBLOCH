function write_WF(x,Psi,Energies,V,Bias,wf_filename,p_filename)

N = length(x);
NrWF = length(Energies);
f = fopen(wf_filename,'w');

en_specifier = '';
wf_specifier = '%.5f ';
for i = 1:NrWF
    en_specifier = [en_specifier '%.5f '];
    wf_specifier = [wf_specifier ' %.5f'];
end
en_specifier = [en_specifier '\n'];
wf_specifier = [wf_specifier ' \n'];


fprintf(f,'%d\n', NrWF);
fprintf(f,'%d\n',N);
fprintf(f,'%f\n',Bias);
fprintf(f,'\n');
fprintf(f,en_specifier,Energies);
fprintf(f,'\n');
fprintf(f,wf_specifier,[x Psi]');
fclose(f);


f = fopen(p_filename,'w');
pot_specifier = '%.5f %.5f \n';
fprintf(f,pot_specifier,[x V]');
fclose(f);


end