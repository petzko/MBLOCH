function [ x_new,err ] = relaxation( M_inv,b,N,x_old)
% Computes one iteration of the relaxation method for the 
% linear system Ax = b, where A has been split into 
% A = M - N ! 
%  the step function will then be : xnew = M^-1*(b+N*x_old)
%   therefore the user should provide M_inv, b, N and the previous value of x !! 
%

    x_new = M_inv*(b+N*x_old);
    err = norm(x_new - x_old);
    
end

