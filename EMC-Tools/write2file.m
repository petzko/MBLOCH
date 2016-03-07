function   write2file( fname, x,y,title ,xname,yname)
%WRITE2FILE Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(fname,'w');
fprintf(fid,'%s \r\n',title); 
fprintf(fid,'########################\r\n');
fprintf(fid,'%s \t,\t %s \r\n',xname,yname);

lx = length(x);
ly = length(y); 

x =reshape(x,lx,1); 
y = reshape(y,ly,1); 

L = min(lx,ly);
for i = 1:L
    fprintf(fid,'%f  ,  %f\r\n',x(i),y(i));
end
    
fclose(fid); 



end

