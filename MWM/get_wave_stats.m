
%
% function [std_w,mu_c,mu2_c,std_c,N]=get_wave_stats(dat,min_pts)
%
% This program calculates Haar wavelet coefficient statistics
% useful for training the MWM.  It uses the DWT subroutine to
% take a Haar wavelet transform.
%
%
% INPUT:  dat     = vector of training data 
%         min_pts = min number of scaling coeffs at coarsest scale
%                   needed for estimating mean + std
% 
% OUTPUT: std_w = vector of stds of wavelet coeffs at different scales
%                 (coarsest scale first)
%         mu_c  = mean of scaling coeff at coarsest level
%	  mu2_c = mean of scaling coeffs^2 
%         std_c = std of scaling coefficient at coarsest level
%         N = number of wavelet scales used
%
% Note:  The program will use as many as data points as possible, but
%        will truncate the last few data point to obtain a feasible 
%        DWT length (ncoarse*2^N).
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

function [std_w,mu_c,mu2_c,std_c,N]=get_wave_stats(dat,min_pts)

N=floor(log2(length(dat)/min_pts));    % # of wavelet scales
if (N<1)
  error('Minimum pts too high for given data length');
end

ncoarse=floor(length(dat)/2^N);        % # of coarse coeffs (>= min_pts)

dat=dat(1:(ncoarse*2^N));              % Discard extra data for DWT
h= [ 1 1]/sqrt(2); g= [-1 1]/sqrt(2);  % Haar wavelet filters
w=haardwt(dat,N);                        % Haar wavelet transform



std_w=zeros(N,1);
for scale=N:-1:1
  first=ncoarse*2^(scale-1)+1;
  last =ncoarse*2^(scale);
  std_w(scale)=sqrt(mean(w(first:last).^2));
end

mu_c=mean(w(1:ncoarse));
std_c=std(w(1:ncoarse));
mu2_c=mean(w(1:ncoarse).^2);



