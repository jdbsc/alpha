
%
% function [sig] = gen_pm_mwm(c,r,mu_csig,std_csig,ncoarse) 
%
% Generates a point-mass MWM trace from wavelet coefficient stds
% and coarest-scale scaling coefficient means and stds.  Note that
% the point-mass MWM may look somewhat artificial, since the data
% will take a finite # of values.
%
% INPUTS: 
%          c        = N_s x 1 vector of location parameters   
%          r        = N_s x 1 vector of weight parameters   
%           
%          ncoarse  = # of coarse-scale coeffs
%                    (Note that choosing ncoarse = M is equivalent
%                     to running this program M times with ncoarse = 1
%                     and concatenating the traces.)
%          mu_csig  = mean of coarse-scale aggregated data
%          std_csig = std of coarse-scale aggregated data
%          (Note that mu_csig and std_csig will be off from the mu_c and
%           std_c as measured by gen_bmwm by a factor of 2^(N_s/2).) 
%
% Outputs:
%          sig     =  (2^N_s * ncoarse) x 1 vector containing point-mass MWM signal.
%
%
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

function [sig] = gen_pm_mwm(c,r,mu_csig,std_csig,ncoarse)

% Keep changing random # seeds.  Remove if want same seed each time.
rand('state',sum(100*clock));   
randn('state',sum(200*clock));

% Find # of scales in Haar wavelet transform
N_s=length(c);

% Synthesize coarsest-scale signal 
sig=(mu_csig+std_csig*randn(1,ncoarse));

if (any(sig<0))
   disp(' ');
   disp('Oops, signal went negative.')
   disp('Normal model for scaling coeffs not appropriate.');
   disp('Maybe (1) try log-normal or some other +ve distribution.')
   disp('      (2) fit data using more wavelet scales.')
   disp(' ');
   sig=[];
   return;
end

for scale=1:N_s
    u=rand(1,ncoarse*2^(scale-1));
    a = c(scale)*(u<r(scale)) - c(scale)*(u>=r(scale)).*(u<2*r(scale)) + 0*(u>=2*r(scale)); 
    sig=[sig.*(1+a)/2;sig.*(1-a)/2];
    sig=reshape(sig,1,ncoarse*2^(scale));
end
  
  
sig=sig';                           % Convert sig to a column vector













