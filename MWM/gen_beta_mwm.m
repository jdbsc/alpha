
%
% function [sig] = gen_beta_mwm(p,mu_c,std_c,ncoarse) 
%
% Generates a betaMWM trace from wavelet coefficient stds
% and coarest-scale scaling coefficient means and stds.
%
%
% INPUTS: 
%          p       = N_s x 1 vector of beta parameters
%                    (N_s is the # of wavelet scales)
%          mu_c    = Mean of coarse-scale scaling coefficients
%          std_c   = Std of coarse-scale scaling coefficients
%          ncoarse = # of coarse-scale scaling coeffs
%                    (Note that choosing ncoarse = M is equivalent
%                     to running this program M times with ncoarse = 1
%                     and concatenating the traces.)
%
% Outputs:
%          sig     =  (2^N_s * ncoarse) x 1 vector containing betaMWM signal.

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

function [sig] = gen_beta_mwm(p,mu_c,std_c,ncoarse)

% Keep changing random # seeds.  Remove if want same seed each time.
rand('state',sum(100*clock));   
randn('state',sum(200*clock));

% Find # of scales in Haar wavelet transform
N_s=length(p);

% Synthesize coarsest-scale signal  
sig=(mu_c+std_c*randn(1,ncoarse))/2^(N_s/2);%产生U00


if (any(sig<0))
   disp(' ');
   disp('Oops, signal went negative.')
   disp('Normal model for scaling coeffs not appropriate.');
   disp('Maybe (1) try log-normal or some other +ve distribution.')
   disp('      (2) fit data using more wavelet scales.')
   disp(' ');
   sig=[];

   % Formulas for a log-normal fit of mean and std
   % lstd=sqrt(log(1+std_c^2/mu_c^2));
   % lmu= log(mu_c)-1/2*lstd^2;
   % sig=exp(lmu+lstd*randn(1,ncoarse))/2^(N_s/2);;

   return;
end

for scale=1:N_s
    a=2*betarnd2(p(scale),p(scale),1,ncoarse*2^(scale-1)) - 1;%产生Ajk
    sig=[sig.*(1+a);sig.*(1-a)];                              %产生Ujk
    sig=reshape(sig,1,ncoarse*2^(scale));
end
  
  
sig=sig';                           % Convert sig to a column vector













