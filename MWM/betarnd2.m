% function [x]=betarnd2(p,q,l,m);
%
% INPUTS: p,q = beta distribution paramters
%         l,m = size of matrix 
% OUTPUT:   x = l x m matrix with random numbers from the specified beta
%               distribution
%
% Beta random number generator based on Chen's BA algorithm
% REF: Continuous univariate distributions (vol 2): Johnson, Kotz, Balakrishnan
%

%
% Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 1999 Rice University
% 
% Permission is granted for use and non-profit distribution providing that this
% notice be clearly maintained. The right to distribute any portion for profit
% or as part of any commercial product is specifically reserved for the author.
%


function [x]=betarnd2(p,q,l,m);

if (nargin<=2)
  l=1;
m=1;
end;

N=l*m;
flag=0;

x=[];
  
   alp=p+q;
  if(min(p,q)<=1) 
    
    B=max(1/p,1/q);
    
  else 
    B=sqrt((alp-2)/(2*p*q-alp));
  end

    G=p+1/B;
    n=N;

    while (flag~=1)    
     u=rand(n,2);
     v=B*log(u(:,1)./(1-u(:,1)));
     w=p*exp(v);
     indices=find((alp*log(alp./(q+w))+G*v-1.3862944>=log(u(:,1).^2.*u(:,2))));
     x=[x (w(indices)./(q+w(indices))).']; 
    if (length(x)==N)
      flag=1;
    else    
      n=N-length(x);
    end  
end


x=reshape(x(:),l,m);






