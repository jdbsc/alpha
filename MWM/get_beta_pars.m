
%
% function [p] = get_beta_pars(std_w,mu_c,std_c,mu2_c) 
%
% Computes betaMWM "p" parameters from wavelet coefficient stds
% and coarest-scale scaling coefficient means and stds.
%
%
% INPUTS: 
%         std_w   = a vector of std's of wavelet coeffs at different scales不同尺度下小波系数的std矢量
%                   (starting from coarsest = 1)
%         mu_c    = mean of scaling coeffs at coarsest scale 最粗标度的标度系数平均值
%         std_c   = std of scaling coeffs at coarsest scale最粗标度的标度系数标准差
%         mu2_c    = mean of scaling coeffs^2 
%
% Outputs:
%         p       = N_s x 1 vector of beta parameters, one for each scale
%                   (e.g. multipliers at scale k are beta(p(k),p(k))
%                         and p(1) corresponds to coarsest scale)
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

function [p] = get_beta_pars(std_w,mu_c,std_c,mu2_c)


% Find # of scales in Haar wavelet transform
N_s=length(std_w);

% Find beta parameter p for coarsest wavelet scale 
p(1)=0.5*( mu2_c/(std_w(1)^2) -1 ); 
if (p<=0) 
  error('Beta parameter p negative => MWM cannot fit data at coarsest scale.');
end

% Compute wavelet energy ratios for calculating p's
eta=(std_w(1:N_s-1)./std_w(2:N_s)).^2;

for scale=1:N_s-1
    p(scale+1)=eta(scale)*.5*(p(scale)+1)-0.5;      % Iterative formula for p's
    if (p<=0) 
       disp(['Scale is ',num2str(scale)]);
       error('Beta parameter p negative => MWM cannot fit data.');
    end
end

  






