
%
% function [p,mu_c,std_c] = train_beta_mwm(sig,min_pts) 
%
% Calculates betaMWM parameters for a given signal.
%
%
%
% INPUTS: 
%         sig     = Column vector containing signal
%         min_pts = min number of scaling coeffs at coarsest scale
%                   needed for estimating mean + std
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

function [p,mu_c,std_c] = train_beta_mwm(sig,min_pts) 

% Calculate wavelet coefficient 2nd-order stats scale-by-scale
[std_w,mu_c,mu2_c,std_c,N]=get_wave_stats(sig,min_pts);

% Convert wavelet coeff stds into p parameters
[p] = get_beta_pars(std_w,mu_c,std_c,mu2_c);




