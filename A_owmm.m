%global y  yy
%global nn
%global nnn
%global epsi T
%%
%导入数据
%[filename1,filepath1]=uigetfile({'*.txt','All Files'},'Select Data File');
%fullpathname=strcat(filepath1,filename1);
fid = fopen('C:\Users\Dell\MATLAB\Projects\Network-Link-Dimensioning-master\loc2-3.txt');
   
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
%figure
c=0:T:nn(end);
% stairs(c,y/T) , xlabel('time (sec)'), ylabel('data rate (bps)'),title(num2str(T)) 
%% 将输入序列y轴变成包到达间隔
for i=1:(length(nn)-1)
    bb(i)=nn(i+1)-nn(i);
  
end

%% 估算数据alpaha分布的参数
y_T=yy';
% estimate parameters
ppp = stblfit(bb','ecf',statset('Display','iter'));
% plot data with fit parameters
% xmax = length(nn);
% H = figure(2);
% set(H,'Position', [517 426 939 410]);
% 
% title('Stable fit to sums of Pareto random variables');
% subplot(1,2,1)
% hold on
% stem(y_T(y_T < xmax),stblpdf(y_T(y_T<xmax),ppp(1),ppp(2),ppp(3),ppp(4),'quick'));
% x = 0:.1:xmax;
% plot(x,stblpdf(x,ppp(1),ppp(2),ppp(3),ppp(4),'quick'),'r-')
% hold off
% xlabel(['\alpha_0 = ',num2str(ppp(1)),'  \beta_0 = ',num2str(ppp(2)),'  \gamma_0 = ',num2str(ppp(3)),'  \delta_0 = ',num2str(ppp(4))]);
% legend('Data','Fit stable density')
% subplot(1,2,2)
% CDF = prctile(y_T,[1:75]);
% cmin = CDF(1);
% cmax = CDF(end);
% x = cmin:.1:cmax;
% estCDF = stblcdf(x,ppp(1),ppp(2),ppp(3),ppp(4));
% plot(CDF,[.01:.01:.75],'b.',x,estCDF,'r-')
% legend('Empirical CDF','Estimated CDF','Location','northwest')% alpha=ppp(1),beta=ppp(2), gammma=ppp(3),delta=ppp(4);


%% 算法主体
%设置迭代层数cen
cen=fix(log2(length(y/T)));
bb=bb';

%生成beta分布的MWm
min_pts=14; ncoarse=16;
[p,mu_c,std_c] = train_beta_mwm(bb,min_pts); 
[MWMgen] = gen_beta_mwm(p,mu_c,std_c,ncoarse);


%生成多层alpha分布的参数，竖向排列
pcan(:,1)=ppp;
[p,mu_c, std_c, N] = train_alpha_mwm(bb,min_pts); 

% x=-1:.001:3;
% plot(x,stblpdf(x,2,0,0.11,0.5))
% hold on
% plot(x,wblpdf(x,0.54,3.6))
% hold off

for i=1:4
    pcan(:,i)=[2;0;0.108;0.49];%2;0;0.11;0.5
  
end
for i=5:8
    pcan(:,i)=[2;0;0.098;0.45];%2;0;0.11;0.5
   
end
for i=9:12
    pcan(:,i)=[2;0;0.089;0.42];%2;0;0.11;0.5
   
end
for i=13:N
    pcan(:,i)=[2;0;0.084;0.41];%2;0;0.11;0.5
   
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
    sd=2*stblrnd(pcan(1,scale),pcan(2,scale),pcan(3,scale),pcan(4,scale),1,ncoarse*2^(scale-1)) - 1; %产生Ajk
    sig=[sig.*(1+sd);sig.*(1-sd)];                                        %产生Ujk
    sig=reshape(sig,1,ncoarse*2^(scale));
end
  
sig=sig'; 

% for si=1:length(sig)
% if sig(si)<0
% 
%     sig(si)=-sig(si);
% end
% 
% end

max(sd)
min(sd)




figure
reciwm=iwm([1:1:length(bb)],bb,16);
figure
plot(bb)
title("original")

figure
h = plot(MWMgen,"DisplayName","testsig1");
ylabel('')
title('')
legend('MWM生成的流量')

figure
h=qqplot(bb,reciwm);
xlabel({'实际分位数','(a) IWM'})
ylabel('理论分位数')
ylim([-0.02 0.12])
legend('实际流量','','IWM')
%title('(a) IWM')
h(3).LineWidth=1.5;

figure
h=qqplot(bb,MWMgen);
legend('实际流量','','MWM')
xlabel({'实际分位数','(b) MWM'})
ylabel('理论分位数')
ylim([0 0.15])
h(3).LineWidth=1.5;

[wowmqq] = wowmqqplot('C:\Users\Dell\MATLAB\Projects\Network-Link-Dimensioning-master\loc2-3.txt');

figure
h=qqplot(bb,wowmqq);
xlabel({'实际分位数','(c) W-OWM'})
ylabel('理论分位数')
legend('实际流量','','W-OWM')
ylim([0 0.15])
h(3).LineWidth=1.5;

figure
h=qqplot(bb,sig);
xlabel({'实际分位数','(d) A-OWM'})
ylabel('理论分位数')
legend('实际流量','','A-OWM')
ylim([0 0.15])
h(3).LineWidth=1.5;

%%
figure
x_pdf=0:1*10^(-4):0.06;
% pd=fitdist(bb','stable');
y_cdf = cdf(pd,x_pdf);
y_qujian=hist(bb,0:1*10^(-3):0.06)/length(bb);
plot(x_pdf,y_cdf*100);
hold on
bar(0:1*10^(-3):0.06,y_qujian*100);
xlabel('数据包到达时间（s）')
ylabel('概率')
ytickformat('percentage')
legend('累计概率','区间概率')
set(gca,'xticklabel',0:0.4:3)


for x_bb=1:1000000-1

bb_d(x_bb)=bb(x_bb);
reciwm_d(x_bb)=reciwm(x_bb);
MWMgen_d(x_bb)=MWMgen(x_bb);
sig_d(x_bb)=sig(x_bb);
wowmg_d(x_bb)=wowmg(x_bb);
end
[dh1,h1,cp1,tauq1] = dwtleader(bb_d);
[dh2,h2,cp2,tauq2] = dwtleader(reciwm_d);
[dh3,h3,cp3,tauq3] = dwtleader(MWMgen_d);
[dh4,h4,cp4,tauq4] = dwtleader(sig_d);
[dh5,h5,cp5,tauq5] = dwtleader(wowmg_d);
figure;
hp = plot(h1,dh1,'-o',h2,dh2,'-^',h3,dh3,'-*',h4,dh4,'k-+',h5,dh5,'m-diamond');
hp(1).MarkerFaceColor = 'b';
hp(2).MarkerFaceColor = 'r';
hp(3).MarkerFaceColor = "#EDB120";
hp(4).MarkerFaceColor = 'k';
hp(5).MarkerFaceColor = 'm';

grid on;
xlabel('h'); ylabel('f(h)');
legend('实际流量','A-OWM','MWM','IWM','W-OWM','Location','NorthEast');
title('多重分形谱');
ylim([-inf 1.2])


%画H值图
% x_bb=1:1:length(bb);
% x_reciwm=1:1:length(reciwm);
% x_MWMgen=1:1:length(MWMgen);
% x_sig=1:1:length(sig);
% 
% 
% rstest(x_bb,bb,1);
% legend('标签流量')
% ax1=gca;
% rstest(x_reciwm,reciwm,1);
% legend('IWM')
% ax2=gca;
% rstest(x_sig,sig,1);
% legend('MWM')
% ax3=gca;
% rstest(x_MWMgen,MWMgen,1);
% legend('A-OWM')
% ax4=gca;
% 
% fnew = figure;
% ax1_copy = copyobj(ax1,fnew);
% ax2_copy = copyobj(ax2,fnew);
% ax3_copy = copyobj(ax3,fnew);
% ax4_copy = copyobj(ax4,fnew);
% subplot(2,2,1,ax1_copy)
% legend('标签流量')
% subplot(2,2,2,ax2_copy)
% legend('IWM')
% subplot(2,2,3,ax3_copy)
% legend('A-OWM')
% subplot(2,2,4,ax4_copy)
% legend('MWM')


% bb_ds=sort(bb_d);
% MWMgen_ds=sort(MWMgen_d);
% sig_ds=sort(sig_d);
% reciwm_ds=sort(reciwm_d);
% R = corrcoef(bb_ds,MWMgen_ds)
% R = corrcoef(bb_d,sig_ds)
% R = corrcoef(bb_d,reciwm_ds)

%% 自相关系数
ee=20;
[acf1,lags1]=autocorr(bb,ee);
[acf2,lags2]=autocorr(reciwm,ee);
[acf3,lags3]=autocorr(MWMgen,ee);
[acf4,lags4]=autocorr(sig,ee);
[acf5,lags5]=autocorr(wowmg,ee);

figure
hp = plot(lags1,acf1,'-o',lags2,acf2,'-^',lags3,acf3,'-*',lags4,acf4,'-+',lags5,acf5,'m-diamond');
hp(1).MarkerFaceColor = 'b';
hp(2).MarkerFaceColor = 'r';
hp(3).MarkerFaceColor = '#77AC30';
hp(4).MarkerFaceColor = 'k';
hp(5).MarkerFaceColor = 'm';
xlabel('延时')
ylabel('自相关系数')
legend('实际流量','MWM','A-OWM','IWM','W-OWM','Location','NorthEast');
