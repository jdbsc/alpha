%
% function [std_w,mu_c,std_c]=get_ss_wave_stats(H,mu_sig,std_sig,N_s);
%
% Generate betaMWM parameters for non-negative signal which is 
% approximately 2nd-order self-similar with Hurst parameter H.
% The parameters are derived from analytic calculations.
% 
%
% INPUTS: H       = Hurst Parameter
%         mu_sig  = mean of process
%         std_sig = std of process
%         N_s     = number of wavelet scales (data length is 2^N_s) 
%
% OUTPUT: std_w = std of wavelet coeffs at different scales
%         (coarsest scale = 1)
%         mu_c  = mean of coarsest-scale scaling coeff
%         std_c = std of coarsest-scale scaling coeff
%
%
%
% Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 1999 Rice University
% 
% Permission is granted for use and non-profit distribution providing that this
% notice be clearly maintained. The right to distribute any portion for profit
% or as part of any commercial product is specifically reserved for the author.
%


function [std_w,mu_c,std_c]=get_ss_wave_stats(H,mu_sig,std_sig,N_s);

std_w=std_sig*sqrt(2-2^(2*H-1));

for count=1:N_s-1
  std_w=[std_w(1)*2^(H-0.5) std_w];
end

mu_c  = mu_sig*2^(N_s/2);
std_c = std_sig*2^(N_s*(H-0.5));





