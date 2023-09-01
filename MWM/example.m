
%the file interarrival.mat contains test data 
load interarrival.mat 
x=interarrival;
 

% Example how to fit data with betaMWM and then synthesize
% Note: first need to set
% x= your_data;
min_pts=15; ncoarse=15;
[p,mu_c,std_c] = train_beta_mwm(x,min_pts); 
[testsig1] = gen_beta_mwm(p,mu_c,std_c,ncoarse);


% Note: first need to set
% x= your_data;
min_pts=10; ncoarse=2;
[c,r,mu_csig,std_csig]=train_pm_mwm(x,min_pts);
[testsig2] = gen_pm_mwm(c,r,mu_csig,std_csig,ncoarse) ;

% Example how to fit data with hybrid beta/point-mass MWM and then synthesize
% Note: first need to set
% x= your_data;
min_pts=15; ncoarse=15;
%x=interarrival(1:15*2^16);
[p,mu_c,std_c,c,r,switchnum] = train_hybrid_mwm(x,min_pts); 
[testsig3] = gen_hybrid_mwm(p,mu_c,std_c,c,r,switchnum,ncoarse);


% Example of using BMWM for synthesizing 2nd-order self-similar traffic
% (traffic with approximately same covariance as fGn with Hurst parameter H)
%
H=0.43; mu_sig=.2; std_sig=.3; N_s=10;ncoarse=15;
[p,mu_c,std_c]=train_ss_beta_mwm(H,mu_sig,std_sig,N_s); 
[testsig4]=gen_beta_mwm(p,mu_c,std_c,ncoarse);

% Example of using the MSQ formula for the tail queue probability for the
% MWM

H=0.8; mu_sig=.2; std_sig=.2; N_s=16;link_speed=.5; queue_sizes=mu_sig*(2.^(1:7)-1);
[p,mu_c,std_c]=train_ss_beta_mwm(H,mu_sig,std_sig,N_s); 
[msq,cdtsq]=mwm_msq(p,mu_c,std_c,link_speed,queue_sizes);









