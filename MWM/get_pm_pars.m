
% function [r,c]=get_pm_pars(mm1,m2)
%
%
% Fits a point-mass model with Pr[A=c]=Pr[A=-c]= r
%                          and Pr[A=0]=1-2r 
% based on desired negative first and second-order moments.
%
%                                                    1+A
% Note the fit is based on matching the moments of  -----
%                                                     2
% to the desired target moments.
%
%
%
% Inputs:
%      mm1 --- A column vector of target -1st moments
%       m2 --- A column vector of target  2nd moments
%
% Outputs:
%        c --- A column vector of pm locations
%        r --- A column vector of pm weights
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


function [r,c]=get_pm_pars(mm1,m2)

% Basic moment-matching formulas
r = - (2-8*m2+4*m2.*mm1-mm1)./(-2*mm1+16*m2);
c =  sqrt((-2+mm1).*(-2+4*r+mm1))./(-2+4*r+mm1); 

% Now check for unacceptable r,c values.  If can't fit both moments,
% match 2nd-moment and make -1st moments as small as possible.
% (Note that we can always make the -1st moment arbitrarily large, 
%  so failure to match implies that the -1st moment of the model is 
%  too large.)

ind = (find ( (r<0)|(r>1)|(imag(c)~=0)|(c>1)|(c<-1)|isinf(r)|isnan(r)));

if(isempty(ind) == 0)
   
   % Conjecture: to minimize -1st moment given a fixed 2nd moment
   %
   %                      (1+c)
   % choose r=0.5 so that  ---- is as far away from 0 as possible.
   %                        2
   % (Not proven, but confident this is true, especially for c close
   % to 0.)

   r(ind)=.5;
   c(ind)=sqrt(4*m2(ind)-1);

   % Check to see that we at least match 2nd moment
   ind2 = ( find((c<-1) | (c>1)) );

   if(isempty(ind2)==0)
      disp(' ')
      disp('Error cannot even match 2nd moments.')
      disp('Matching as closely as possible.')
      disp(' ')
      r(ind(ind2))=0.5*ones(size(ind2));
      c(ind(ind2))=ones(size(ind2));
   end
end





