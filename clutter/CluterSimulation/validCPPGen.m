%% 
% Rayleigh
clear all;
ydata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\rayleigh.txt");
fr=12e3;
sigmav=40;
sigmaf=2*sigmav*6e9/3e8;
sigmac=1.2 ;           %杂波的标准差
fftNum=2^ceil(log2(length(ydata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure;
subplot(211);plot(real(ydata)); %瑞利分布杂波实部
title('瑞利杂波时域波形，实部');  
subplot(212);plot(imag(ydata)); %瑞利分布杂波虚部
title('瑞利杂波时域波形，虚部')


num=100;                   %求概率密度函数的参数
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);   
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %用直方图估计的概率密度函数
xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %概率密度函数理论值
figure;
plot(xaxisl,xpdf1);               %做出仿真结果的概率密度函数曲线
hold on;plot(xaxisl,th_val,'r:'); %做出理论概率密度函数曲线
title('c++生成 瑞利 杂波幅度分布');
xlabel('幅度');ylabel('概率密度');

signal=ydata;
signal=signal-mean(ydata);
figure;
M=256;
psd_dat=fftshift(fft(signal, fftNum));
freqx=f;
plot(freqx,abs(psd_dat));
title('杂波频谱'),xlabel('频率/Hz'),
ylabel('功率谱密度');
%%%%%%%% 理想高斯曲线%%%%%%%%%%%%
powerf=exp((-abs(f)/sigmaf))*200;
y=freqx*3e8/2/6e9;
hold on;
plot(freqx,powerf,':r');
%% 
% Lognormal
clear all;
xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\lognormal.txt");

fr=12e3;
sigmav=40;     %杂波方差
sigmaf=2*sigmav*6e9/3e8;
muc=1.5;
sigmac=0.6;
fftNum=2^ceil(log2(length(xdata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure,plot(abs(xdata));xlabel('时间');ylabel('幅度');title('对数正态分布杂波时域波形');
num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=lognpdf(xaxis1,log(muc),sigmac);

figure;plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r');
title('c++生成 对数正态 杂波幅度分布');xlabel('幅度');ylabel('概率密度');

signal=xdata;
signal=signal-mean(signal);

figure,M=128;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,psd_dat);title('杂波频谱');xlabel('频率(Hz)');ylabel('功率谱密度');
powerf=exp(-freqx.^2/(2*sigmaf.^2))*150;
hold;plot(freqx,abs(powerf),':r');
%% 
%weibull
clear all;
xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\weibull.txt");

fr=12e3;
sigmav=40;     %杂波方差
sigmaf=2*sigmav*6e9/3e8;
p=1.5;%形状参数
q=2.2;%尺度参数
fftNum=2^ceil(log2(length(xdata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure,plot(abs(xdata));xlabel('时间');ylabel('幅度');title('韦布尔分布杂波时域波形');

num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=p*(xaxis1.^(p-1)).*exp(-(xaxis1/q).^p)./(q.^p);

figure;plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r');
title('c++生成 韦布尔 杂波幅度分布');xlabel('幅度');ylabel('概率密度');

signal=xdata;
signal=signal-mean(signal);

figure,M=256;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,abs(psd_dat));title('杂波频谱');xlabel('频率(Hz)');ylabel('功率谱密度');
powerf=exp(-freqx.^2/(2*sigmaf.^2))*200;
hold;plot(freqx,powerf,':r');
%%
%k分布
clear all;
clc;
ydata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\k.txt");
xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\kxdata.txt");

fr=12e3;
sigmav=40;     %杂波方差
sigmaf=2*sigmav*6e9/3e8;
vmuc=2;
fftNum=2^ceil(log2(length(xdata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure,subplot(2,1,1),plot(real(ydata)), title('K分布杂波时域波形，实部');
subplot(2,1,2),plot(imag(ydata)),title('K分布杂波时域波形，虚部');
%%%%%%%%求概率密度函数的参数%%%%%%%%%%%%%%%%%%
num=100;
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));
alpha=sqrt(std(xdata).^2./(2*vmuc));%std()算xdata标准差

% xpdf1=getnpdf(abs(xdata),num,maxdat,mindat);
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=2*((xaxis1/(2*alpha)).^vmuc).*besselk((vmuc-1),xaxis1/alpha)./(alpha*gamma(vmuc));
figure,plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r'),title('c++ k分布 杂波幅度分布'),xlabel('幅度'),
ylabel('概率密度');
signal=ydata;
signal=signal-mean(ydata);
M=256;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

figure,plot(freqx,abs(psd_dat));
title('杂波频谱'),xlabel('频率/Hz'),
ylabel('功率谱密度');
%%%%%%%% 理想高斯曲线%%%%%%%%%%%%
powerf=exp(-freqx.^2/(2*sigmaf.^2))*200;
hold;
plot(freqx,powerf,':r');


% clear;
% data = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\save\k.txt");
% xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\save\kdata.txt");
% ydata = data;
% vmuc = 2;
% figure,subplot(2,1,1),plot(real(ydata)), title('K分布杂波时域波形，实部');
% subplot(2,1,2),plot(imag(ydata)),title('K分布杂波时域波形，虚部');
% %%%%%%%%求概率密度函数的参数%%%%%%%%%%%%%%%%%%
% num=100;
% maxdat=max(abs(ydata));
% mindat=min(abs(ydata));
% NN=hist(abs(ydata),num);
% xpdf1=num*NN/((sum(NN))*(maxdat-mindat));
% xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
% alpha=sqrt(std(xdata).^2./(2*vmuc));%std()算xdata标准差
% th_val=lognpdf(xaxis1,xpdf1);
% 
% % xpdf1=getnpdf(abs(xdata),num,maxdat,mindat);
% xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
% th_val=2*((xaxis1/(2*alpha)).^vmuc).*besselk((vmuc-1),xaxis1/alpha)./(alpha*gamma(vmuc));
% figure,plot(xaxis1,xpdf1);
% hold,plot(xaxis1,th_val,':r'),title('c++ k分布 杂波幅度分布'),xlabel('幅度'),
% ylabel('概率密度');