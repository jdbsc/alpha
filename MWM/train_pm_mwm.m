
%
% function [c,r,mu_c,std_c] = train_pm_mwm(sig,min_pts) 
%
% Calculates point-massMWM  parameters for a given signal.
% The point-mass function has a mass of r at +- c and 
% 1-2r at 0.
%
% INPUTS: 
%         sig     = Column vector containing signal
%         min_pts = min number of scaling coeffs at coarsest scale
%                   needed for estimating mean + std
%
%
% Outputs:
%         c       = N_s x 1 vector of point-mass locations, one for each level
%         r       = N_s x 1 vector of point-mass weights
%         mu_c    = mean of data at coarsest level
%         std_c   = std of data at coarsest level
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

function [c,r,mu_csig,std_csig] = train_pm_mwm(sig,min_pts) 

% Calculate moment ratios (2nd and -1st) of data
% (-1st moment is good if you want to fit small values
%  you can also use higher-order moments)
%
[mm1_ratios,m2_ratios,mu_csig,std_csig]=get_moment_ratios(sig,min_pts);


% Convert moment ratios into r,c parameters
[r,c] = get_pm_pars(mm1_ratios,m2_ratios);

