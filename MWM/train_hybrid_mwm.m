
%
% function [p,mu_c,std_c,c,r,switchnum] = train_hybrid_mwm(sig,min_pts) 
%
% Calculates hybrid MWM parameters for a given signal.  The beta is used 
% until error in the negative 1st moment (when compared to actual)
% exceeds an error tolerance set by the variable "error_tol".
%
%
% INPUTS: 
%         sig     = Column vector containing signal
%         min_pts = min number of scaling coeffs at coarsest scale
%                   needed for estimating mean + std
%
% Outputs:
%         p         = N_s x 1 vector of beta parameters, one for each scale
%                   (e.g. multipliers at scale k are beta(p(k),p(k))
%                         and p(1) corresponds to coarsest scale)
%         mu_c      = mean of scaling coefficients at coarsest scale
%         std_c     = std of scaling coefficients at coarsest scale
%         c         = N_s x 1 vector of point-mass locations
%         r         = N_s x 1 vector of point-mass weights
%         switchnum = Flag indicating when to switch from beta to point-mass
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

function [p,mu_c,std_c,c,r,switchnum] = train_hybrid_mwm(sig,min_pts )

error_tol = 0.05;   % error_tolerance for -1st moment


% Fit beta MWM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate wavelet coefficient 2nd-order stats scale-by-scale
[std_w,mu_c,mu2_c,std_c,N]=get_wave_stats(sig,min_pts);
% Convert wavelet coeff stds into p parameters
[p] = get_beta_pars(std_w,mu_c,std_c,mu2_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit point-mass MWM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate moment ratios
[mm1_ratios,m2_ratios,mu_csig,std_csig]=get_moment_ratios(sig,min_pts);
% Convert moment ratios into r,c parameters
[r,c] = get_pm_pars(mm1_ratios,m2_ratios);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure where to switch from beta to point-mass
switchnum=N+1;  
for scale=1:N

  if (p(scale)>1)
      beta_mm1(scale) = (2*p(scale)-1)/(p(scale)-1);  % Theor -1st moment
  else
      beta_mm1(scale) = Inf;
  end

  if ((abs(beta_mm1(scale)-mm1_ratios(scale))/mm1_ratios(scale) > error_tol) )
     switchnum=scale;
     break
  end
end    
     





