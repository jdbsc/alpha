function [sig]=wowmqqplot(filepath)

%fid = fopen('C:\Users\Dell\MATLAB\Projects\Network-Link-Dimensioning-master\loc2-3.txt');
fid = fopen(filepath);   
a = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.

a=a' ;
nn=a(:,1);  % time vector
nnn=a(:,2)*8 ; % packets length 
x=str2num('0.01');              % times of excution   
epsi=str2num('0.1');              % 

%% 将导入数据画图
y=0;
for kk=1:length(x)
    sto=0;y=0;z=1;
    T=x(kk);
    
    for i=1:length(nn)
        if nn(i)>=x(kk)
            sto(z)=i;
            x(kk)=x(kk)+T;
            z=z +1;
        end
    end
    sto=[1 sto length(nn)+1];
    
    for j=1:length(sto)-1  
        y(j)=sum(nnn(sto(j):sto(j+1)-1)); 
    end
    yy=y/T;  %% yy==>bit/sec
end
%figure(1)
%c=0:T:nn(end);
% stairs(c,y/T) , xlabel('time (sec)'), ylabel('data rate (bps)'),title(num2str(T)) 
%% 将输入序列y轴变成包到达间隔
for i=1:(length(nn)-1)
    bb(i)=nn(i+1)-nn(i);
  
end

%% 算法主体
cen=fix(log2(length(y/T)));

min_pts=14; ncoarse=16;

%生成多层weibull分布的参数，竖向排列
pcan=zeros(2,1);

for u=1:length(bb)
    if bb(u)==0
        bb(u)=0.000001;
    end
end

pd=fitdist(bb','weibull');

[p,mu_c, std_c, N] = train_alpha_mwm(bb,min_pts); 
for i=1:4
    pcan(:,i)=[0.538,3.6];%2;0;0.11;0.5
   
end
for i=5:8
    pcan(:,i)=[0.538,3.65];%2;0;0.11;0.5
   
end
for i=9:12
    pcan(:,i)=[0.538,3.65];%2;0;0.11;0.5
   
end
for i=13:N
    pcan(:,i)=[0.538,3.79];%2;0;0.11;0.5
   
end

% Keep changing random # seeds.  Remove if want same seed each time.
rand('state',sum(100*clock));   
randn('state',sum(200*clock));

% Find # of scales in Haar wavelet transform
N_s=length(pcan(1,:));

% Synthesize coarsest-scale signal  
sig=(mu_c+std_c*randn(1,ncoarse))/2^(N_s/2); %产生U0k

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
    sd=2*wblrnd(pcan(1,scale),pcan(2,scale),1,ncoarse*2^(scale-1))-1; %产生Ajk
    sig=[sig.*(1+sd);sig.*(1-sd)];                                        %产生Ujk
    sig=reshape(sig,1,ncoarse*2^(scale));
end
  
sig=sig'; 
% ddddd=wblrnd(pd.A,pd.B,1,1000000);


end
