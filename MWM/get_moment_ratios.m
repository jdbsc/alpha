
% [mm1_ratios,m2_ratios,mu_c,std_c,N]=get_moment_ratios(sig,min_pts)
%
% This routine calculates the moment ratios for data at different
% resolution levels.  These ratios are useful fitting a point-mass
% density for the MWM.  
% 
% To keep things simple, moments are calculated on raw aggregated 
% data (e.g. agg(x)= x(1:2:N-1)+x(2:2:N) ), so there is no use of 
% wavelet transform or the 2^(j/2) type scaling factors needed for 
% scaling coefficients.
%
%
%
% INPUTS: 
%         sig     = Column vector containing signal
%         min_pts = min number of data coeffs at coarsest scale
%                   needed for estimating mean + std
%
% Outputs:
%         mm1_ratios - (nlevels-1)x1 ratios of negative first moments 
%                      mm1_ratios(1) is ratio at coarsest scale               
%         m2_ratios  - (nlevels-1)x1 ratios of second moments
%         mu_csig    - mean of signal at coarsest aggregation level
%         std_csig   - variance of signal coarsest aggregation level
%

function [mm1_ratios,m2_ratios,mu_csig,std_csig]=get_moment_ratios(sig,min_pts)

nlevels=floor(log2(length(sig)/min_pts));  
if (nlevels<1)
  error('Minimum pts too high for given data length');
end

ncoarse=floor(length(sig)/2^nlevels);
sig=sig(1:ncoarse*2^nlevels);   % Truncate data to keep moments self-consistent

m2(nlevels+1)=mean(sig.^2);
mm1(nlevels+1)=mean(1./sig);

for i=nlevels:-1:1
   sig=sig(1:2:length(sig)-1)+sig(2:2:length(sig));  % Aggregate signal

   m2(i)=mean(sig.^2);
   mm1(i)=mean(1./sig);

   m2_ratios(i)=m2(i+1)/m2(i);                 % Calculate moment ratios finer/coarser
   mm1_ratios(i)=mm1(i+1)/mm1(i);
end

 
std_csig=std(sig);   mu_csig=mean(sig); 





