%
% function [p,mu_c,std_c] = train_ss_beta_mwm(H,mu_sig,std_sig,N_s) 
%
% Calculates betaMWM parameters for a second-order self-similar signal.
%
%
% INPUTS: 
%         H       = H Hurst parameter.  Covariance r(n) ~ n^{2-2H}, 0.5 < H < 1
%         mu_sig  = mean of signal
%         std_sig = std of signal
%         N_s     = # of desired wavelet scales
%
% Outputs:
%         p       = N_s x 1 vector of beta parameters, one for each scale
%                   (e.g. multipliers at scale k are beta(p(k),p(k))
%                         and p(1) corresponds to coarsest scale)
%         mu_c    = mean of scaling coefficients at coarsest scale
%         std_c   = std of scaling coefficients at coarsest scale
%
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

function [p,mu_c,std_c] = train_ss_beta_mwm(H,mu_sig,std_sig,N_s) 

% Calculate wavelet coefficient 2nd-order stats scale-by-scale
[std_w,mu_c,std_c]=get_ss_wave_stats(H,mu_sig,std_sig,N_s) ;

mu2_c=mu_c^2+std_c^2;
% Convert wavelet coeff stds into p parameters
[p] = get_beta_pars(std_w,mu_c,std_c,mu2_c);




