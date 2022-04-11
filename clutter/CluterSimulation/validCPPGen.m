%% 
% Rayleigh
clear all;
ydata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\rayleigh.txt");
fr=12e3;
sigmav=40;
sigmaf=2*sigmav*6e9/3e8;
sigmac=1.2 ;           %�Ӳ��ı�׼��
fftNum=2^ceil(log2(length(ydata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure;
subplot(211);plot(real(ydata)); %�����ֲ��Ӳ�ʵ��
title('�����Ӳ�ʱ���Σ�ʵ��');  
subplot(212);plot(imag(ydata)); %�����ֲ��Ӳ��鲿
title('�����Ӳ�ʱ���Σ��鲿')


num=100;                   %������ܶȺ����Ĳ���
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);   
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %��ֱ��ͼ���Ƶĸ����ܶȺ���
xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %�����ܶȺ�������ֵ
figure;
plot(xaxisl,xpdf1);               %�����������ĸ����ܶȺ�������
hold on;plot(xaxisl,th_val,'r:'); %�������۸����ܶȺ�������
title('c++���� ���� �Ӳ����ȷֲ�');
xlabel('����');ylabel('�����ܶ�');

signal=ydata;
signal=signal-mean(ydata);
figure;
M=256;
psd_dat=fftshift(fft(signal, fftNum));
freqx=f;
plot(freqx,abs(psd_dat));
title('�Ӳ�Ƶ��'),xlabel('Ƶ��/Hz'),
ylabel('�������ܶ�');
%%%%%%%% �����˹����%%%%%%%%%%%%
powerf=exp((-abs(f)/sigmaf))*200;
y=freqx*3e8/2/6e9;
hold on;
plot(freqx,powerf,':r');
%% 
% Lognormal
clear all;
xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\lognormal.txt");

fr=12e3;
sigmav=40;     %�Ӳ�����
sigmaf=2*sigmav*6e9/3e8;
muc=1.5;
sigmac=0.6;
fftNum=2^ceil(log2(length(xdata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure,plot(abs(xdata));xlabel('ʱ��');ylabel('����');title('������̬�ֲ��Ӳ�ʱ����');
num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=lognpdf(xaxis1,log(muc),sigmac);

figure;plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r');
title('c++���� ������̬ �Ӳ����ȷֲ�');xlabel('����');ylabel('�����ܶ�');

signal=xdata;
signal=signal-mean(signal);

figure,M=128;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,psd_dat);title('�Ӳ�Ƶ��');xlabel('Ƶ��(Hz)');ylabel('�������ܶ�');
powerf=exp(-freqx.^2/(2*sigmaf.^2))*150;
hold;plot(freqx,abs(powerf),':r');
%% 
%weibull
clear all;
xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\weibull.txt");

fr=12e3;
sigmav=40;     %�Ӳ�����
sigmaf=2*sigmav*6e9/3e8;
p=1.5;%��״����
q=2.2;%�߶Ȳ���
fftNum=2^ceil(log2(length(xdata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure,plot(abs(xdata));xlabel('ʱ��');ylabel('����');title('Τ�����ֲ��Ӳ�ʱ����');

num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=p*(xaxis1.^(p-1)).*exp(-(xaxis1/q).^p)./(q.^p);

figure;plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r');
title('c++���� Τ���� �Ӳ����ȷֲ�');xlabel('����');ylabel('�����ܶ�');

signal=xdata;
signal=signal-mean(signal);

figure,M=256;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,abs(psd_dat));title('�Ӳ�Ƶ��');xlabel('Ƶ��(Hz)');ylabel('�������ܶ�');
powerf=exp(-freqx.^2/(2*sigmaf.^2))*200;
hold;plot(freqx,powerf,':r');
%%
%k�ֲ�
clear all;
clc;
ydata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\k.txt");
xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\output\kxdata.txt");

fr=12e3;
sigmav=40;     %�Ӳ�����
sigmaf=2*sigmav*6e9/3e8;
vmuc=2;
fftNum=2^ceil(log2(length(xdata)));
f=(-fftNum/2:1:fftNum/2-1)/fftNum*fr;

figure,subplot(2,1,1),plot(real(ydata)), title('K�ֲ��Ӳ�ʱ���Σ�ʵ��');
subplot(2,1,2),plot(imag(ydata)),title('K�ֲ��Ӳ�ʱ���Σ��鲿');
%%%%%%%%������ܶȺ����Ĳ���%%%%%%%%%%%%%%%%%%
num=100;
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));
alpha=sqrt(std(xdata).^2./(2*vmuc));%std()��xdata��׼��

% xpdf1=getnpdf(abs(xdata),num,maxdat,mindat);
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=2*((xaxis1/(2*alpha)).^vmuc).*besselk((vmuc-1),xaxis1/alpha)./(alpha*gamma(vmuc));
figure,plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r'),title('c++ k�ֲ� �Ӳ����ȷֲ�'),xlabel('����'),
ylabel('�����ܶ�');
signal=ydata;
signal=signal-mean(ydata);
M=256;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

figure,plot(freqx,abs(psd_dat));
title('�Ӳ�Ƶ��'),xlabel('Ƶ��/Hz'),
ylabel('�������ܶ�');
%%%%%%%% �����˹����%%%%%%%%%%%%
powerf=exp(-freqx.^2/(2*sigmaf.^2))*200;
hold;
plot(freqx,powerf,':r');


% clear;
% data = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\save\k.txt");
% xdata = dlmread("E:\Documents\Projects\Cpp\SignalProcessFunction\save\kdata.txt");
% ydata = data;
% vmuc = 2;
% figure,subplot(2,1,1),plot(real(ydata)), title('K�ֲ��Ӳ�ʱ���Σ�ʵ��');
% subplot(2,1,2),plot(imag(ydata)),title('K�ֲ��Ӳ�ʱ���Σ��鲿');
% %%%%%%%%������ܶȺ����Ĳ���%%%%%%%%%%%%%%%%%%
% num=100;
% maxdat=max(abs(ydata));
% mindat=min(abs(ydata));
% NN=hist(abs(ydata),num);
% xpdf1=num*NN/((sum(NN))*(maxdat-mindat));
% xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
% alpha=sqrt(std(xdata).^2./(2*vmuc));%std()��xdata��׼��
% th_val=lognpdf(xaxis1,xpdf1);
% 
% % xpdf1=getnpdf(abs(xdata),num,maxdat,mindat);
% xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
% th_val=2*((xaxis1/(2*alpha)).^vmuc).*besselk((vmuc-1),xaxis1/alpha)./(alpha*gamma(vmuc));
% figure,plot(xaxis1,xpdf1);
% hold,plot(xaxis1,th_val,':r'),title('c++ k�ֲ� �Ӳ����ȷֲ�'),xlabel('����'),
% ylabel('�����ܶ�');