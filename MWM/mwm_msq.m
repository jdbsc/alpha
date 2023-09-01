% function [msq,cdtsq]=mwm_msq(p,mu_c,std_c,link_speed,queue_sizes)
%
% Computes approximations (the multiscale queuing formula (MSQ) and the
% critical dyadic time scale queue (CDTSQ)) to P(Q>q), the probability of 
% the queue size exceeding "q" for an infinite queue with constant service 
% rate and fed with traffic from the Multifractal Wavelet Model (MWM).
%
% INPUTS: p          = vector with beta parameters of the MWM model with
%                      p(1) belonging to coarsest scale
%         mu_c       = mean of coarsest scale scaling coefficient
%         std_c      = standard deviation of coarsest scale scaling coefficient
%         link_speed = number of bytes emptied from the queue every time
%                      unit corresponding to the finest time scale.
%         queue_sizes= vector with the different values of "q" for which to 
%                      compute approximations for P(Q>q) 
%
% OUTPUT: msq        = vector with multiscale queuing formula approx. to
%                      P(Q>q)
%         cdtsq      = vector with critical dyadic time scale queue approx. to
%                      P(Q>q)
%
% Example: 
%  H=0.8;
%  mu_sig=100;
%  std_sig=50;
%  N_s=16;
%  link_speed=130;
%  queue_sizes=1:10:500;
%  [p,mu_c,std_c] = train_ss_beta_mwm(H,mu_sig,std_sig,N_s)
%  [msq,cdtsq]=mwm_msq(p,mu_c,std_c,link_speed,queue_sizes)
%
%
% Reference: V. Ribeiro, R. Riedi, M. S. Crouse, R. Baraniuk, "Multiscale 
% queuing analysis of long range dependent traffic," Proceeding of IEEE
% INFOCOM 2000, Tel Aviv, Israel.
% URL: www.dsp.rice.edu/publications
%
%
% Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 2001 Rice University
% 
% Permission is granted for use and non-profit distribution providing that this
% notice be clearly maintained. The right to distribute any portion for profit
% or as part of any commercial product is specifically reserved for the author.
%

function [msq,cdtsq]=mwm_msq(p,mu_c,std_c,link_speed,queue_sizes)

p_c=0.5*((mu_c/std_c)^2-1);
q_c=p_c;


N=length(p);% total scales in Q expression
t=1:N+1;
p=[p_c p];

K=2*mu_c*2^(N/2);
% Approximating the product of beta random variables by another beta
% random variable, we model data at different time scales as beta random 
% variables.


for count=1:N,
  S(count)=2^(-count-1);
  T(count)=S(count)*prod((1+p(1:count+1))./(1+2*p(1:count+1)));
end;

a=S.*(S-T).*((T-S.^2).^-1);
a=[p_c a]; %first beta parameters at different time scales.
a=repmat(a,length(queue_sizes),1);
q=(1-S).*(S-T).*((T-S.^2).^-1);
q=[q_c q]; %second beta parameters at different time scales.
q=repmat(q,length(queue_sizes),1);

bounds=(1/K)*(link_speed*2.^(N+1-t));
bounds=repmat(bounds,length(queue_sizes),1);
queue_sizes=repmat(queue_sizes,N+1,1);
queue_sizes=queue_sizes.';
bounds=bounds+queue_sizes/K;
tmp=find(bounds>1);
bounds(tmp)=1; % betainc requires values <= 1.

[m,n]=size(bounds);
for count=1:m,
  for l=1:n,
    P(count,l)=betainc(bounds(count,l),a(count,l),q(count,l));
  end;
end;
P=P.';
%P=prod((betainc(bounds,a,q)).');
cdtsq=1-min(P);
P=prod(P);
msq=1-P;


