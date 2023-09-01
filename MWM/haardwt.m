function [w]=haardwt(x,nlevels)

% 
% Input:     x - Data of length c*2^nlevels, c an arbitrary positive integer.
%            nlevels - Number of levels (scales) in Haar transform
%
% Output:    w - Wavelet coefficients, ordered as in Rice wavelets toolbox dwt.m.
%                
% 
% Example:  x=randn(1024,1);  [w]=haardwt(x,10);
% 
% File Name: haardwt.m
% Last Modification Date: %G%	%U%
% Current Version: %M%	%I%
% File Creation Date: Tue Feb 23 16:44:37 1999
% Author: Matthew Crouse  <mcrouse@ece.rice.edu>
% 
% Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 1999 Rice University
% 
% Permission is granted for use and non-profit distribution providing that this
% notice be clearly maintained. The right to distribute any portion for profit
% or as part of any commercial product is specifically reserved for the author.
% 
% Change History:
% 

w = zeros(size(x));
L=length(x);

npts=L;

for i=1:nlevels
  lowind=round(npts/2+1); highind=round(npts);
  w(lowind:highind)=(x(1:2:npts-1)-x(2:2:npts))/sqrt(2); % HP/wave coeffs 
  x = (x(1:2:npts-1)+x(2:2:npts))/sqrt(2);               % LP/scale coeffs
  npts=round(npts/2);
end
 

w(1:lowind-1)=x;
