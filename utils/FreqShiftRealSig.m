function [ s_shift ] = FreqShiftRealSig(  sig,dt,skipctr,dv, times)
%SHIFTREALSIG Summary of this function goes here
%   Detailed explanation goes here


    N = length(sig);
    sig = reshape(sig,[N 1]);
% 01) take the complex analytic of E
    s_analytic = hilbert(sig);
% 02) shift the complex analytiic with 2*pi*dv 
    s_analytic = s_analytic.*exp(1i*2*pi*dv*times);
% 03) Take the real part of the shifted product... 
    s_shift = real(s_analytic);
    
    
end

