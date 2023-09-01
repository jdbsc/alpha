
%
% function [sig] = gen_hybrud_mwm(p,mu_c,std_c,c,r,switchnum,ncoarse) 
%
% Generates a hybrid MWM trace (uses beta in coarse scales and point-mass 
% in fine scales). 
%
% Before calling this routine, you must set the beta parameters 
% via "train_hybrid_bmwm".  See that routine for more details.
%
% INPUTS: 
%          p         = N_s x 1 vector of beta parameters   
%                     (N_s is the # of wavelet scales)
%          mu_c      = Mean of coarse-scale scaling coefficients
%          std_c     = Std of coarse-scale scaling coefficients
%          c         = N_s x 1 vector of point-mass locations
%          r         = N_s x 1 vector of point-mass weights
%          switchnum = Flag indicating when to switch from beta to point-mass
%          ncoarse   = # of coarse-scale scaling coeffs
%                    (Note that choosing ncoarse = M is equivalent
%                     to running this program M times with ncoarse = 1
%                     and concatenating the traces.)
%
% Outputs:
%          sig     =  (2^N_s * ncoarse) x 1 vector containing hybrid MWM signal.
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

function [sig] = gen_hybrid_mwm(p,mu_c,std_c,c,r,switchnum,ncoarse)

% Keep changing random # seeds.  Remove if want same seed each time.
rand('state',sum(100*clock));   
randn('state',sum(200*clock));

% Find # of scales in Haar wavelet transform
N_s=length(p);

% Synthesize coarsest-scale signal  
sig=(mu_c+std_c*randn(1,ncoarse))/2^(N_s/2);

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

    % Decide whether to use beta or point-mass
    if (scale<switchnum)
        a=2*betarnd2(p(scale),p(scale),1,ncoarse*2^(scale-1)) - 1;        
    else
        u=rand(1,ncoarse*2^(scale-1));
        a= c(scale)*(u<r(scale))- c(scale)*(u>=r(scale)).*(u<2*r(scale)) + 0*(u>=2*r(scale)); 
    end

    sig=[sig.*(1+a);sig.*(1-a)];
    sig=reshape(sig,1,ncoarse*2^(scale));
end
  
  
sig=sig';                           % Convert sig to a column vector













